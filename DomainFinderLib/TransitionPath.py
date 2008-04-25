from MMTK import *
from MMTK.PDB import PDBConfiguration
from MMTK.Proteins import Protein
from MMTK.ForceFields import CalphaForceField
from MMTK.NormalModes import NormalModes, SparseMatrixSubspaceNormalModes
from MMTK.FourierBasis import FourierBasis, countBasisVectors
from MMTK.Random import uniform, gaussian, randomParticleVector
from MMTK.Trajectory import Trajectory, SnapshotGenerator, TrajectoryOutput
from Scientific import N
import os, string
from random import randrange

def massWeightedNormalVector(v):
    return v/N.sqrt(v.massWeightedDotProduct(v))

# An optimized universe that eliminates the need to recalculate
# the boundary box for each normal mode calculation.
class TransitionPathUniverse(InfiniteUniverse):

    def __init__(self, *args, **kwargs):
        apply(InfiniteUniverse.__init__, (self,)+args, kwargs)
        self.cached_bbox = None
              
    def boundingBox(self, conf = None):
        if self.cached_bbox is None:
            self.setBoundingBox(None, None)
        return self.cached_bbox

    def setBoundingBox(self, conf1, conf2):
        min1, max1 = InfiniteUniverse.boundingBox(self, conf1)
        min2, max2 = InfiniteUniverse.boundingBox(self, conf2)
        min = Vector(N.minimum(min1.array, min2.array))
        max = Vector(N.maximum(max1.array, max2.array))
        self.cached_bbox = min, max

    def description(self, objects=None, index_map = None):
        d = InfiniteUniverse.description(objects, index_map)
        return string.replace(d, self.__class__.__name__,
                              InfiniteUniverse.__name__, 1)

class TransitionPathStep:

    def __init__(self, conf, nmodes):
        self.conf = copy(conf)
        self.universe = conf.universe
        self.np = self.universe.numberOfAtoms()
        self.nmodes = nmodes
        self.calcModes()

    def distance(self, other):
        return (other.conf-self.conf).norm()/N.sqrt(self.np)

    def differenceVectorTo(self, other):
        return other.conf - self.conf

    def moveBy(self, direction, distance):
        vector = direction.scaledToNorm(distance*N.sqrt(self.np))
        self.universe.setConfiguration(self.conf + vector)
        self.fixStructure()
        return TransitionPathStep(self.universe.configuration(), self.nmodes)

    def modeWeighting(self, vector):
        new = ParticleVector(self.universe)
        for n in range(6, len(self.modes)/2):
            m = self.modes[n]
            weight = 1./m.frequency
            new = new + weight*m.massWeightedDotProduct(vector)*m
        return new

    def stepPenalties(self, other):
        d = massWeightedNormalVector(self.differenceVectorTo(other))
        global_term = 0.
        for n in range(6):
            p = self.modes[n].massWeightedDotProduct(d)
            global_term = global_term + p**2
        deformation_term = 0.
        for n in range(6, len(self.modes)):
            m = self.modes[n]
            p = m.frequency*m.massWeightedDotProduct(d)
            deformation_term = deformation_term + p**2
        return global_term, deformation_term

    def calcModes(self):
        self.universe.setConfiguration(self.conf)
        natoms = self.universe.numberOfAtoms()
        nmodes = min(self.nmodes, 3*natoms)
        if nmodes > natoms:
            nmodes = 3*natoms
        if nmodes == 3*natoms:
            modes = NormalModes(self.universe, None)
            modes.cutoff = None
        else:
            p1, p2 = self.universe.boundingBox()
            cutoff_max = (p2-p1).length()
            cutoff = 0.5*cutoff_max
            nmodes_opt = nmodes
            nmodes = countBasisVectors(self.universe, cutoff)
            while nmodes > nmodes_opt:
                cutoff = cutoff + 0.1
                if cutoff > cutoff_max:
                    cutoff = cutoff_max
                    break
                nmodes = countBasisVectors(self.universe, cutoff)
            while nmodes < nmodes_opt:
                cutoff = cutoff - 0.1
                if cutoff < 0.1:
                    cutoff = 0.1
                    break
                nmodes = countBasisVectors(self.universe, cutoff)
            basis = FourierBasis(self.universe, cutoff)
            basis.may_modify = 1
            modes = SparseMatrixSubspaceNormalModes(self.universe, basis, None)
            modes.cutoff = cutoff
        self.modes = modes

    def fixStructure(self):
        conf = self.universe.configuration()
        protein = self.universe[0]
        for chain in protein:
            for i in range(1, len(chain)):
                d = chain[i].neighbour_distance
                r1 = conf[chain[i-1].peptide.C_alpha]
                r2 = conf[chain[i].peptide.C_alpha]
                conf[chain[i].peptide.C_alpha] = r1+d*(r2-r1).normal()


class TransitionPath:

    def __init__(self, universe, conf1, conf2, step_length, nmodes):
        self.universe = universe
        self.conf1 = conf1
        self.conf2 = conf2
        self.step_length = step_length
        self.nmodes = nmodes
        self.setDistances()
        self.path = [[TransitionPathStep(self.conf1, self.nmodes)],
                     [TransitionPathStep(self.conf2, self.nmodes)]]

    def setDistances(self):
        protein = self.universe[0]
        for chain in protein:
            for i in range(1, len(chain)):
                d = (chain[i].position()-chain[i-1].position()).length()
                chain[i].neighbour_distance = d

    def view(self):
        viewSequence(self.universe, map(lambda s: s.conf, self.path))

    def writeBestToTrajectory(self, filename, comment):
        self._writeToTrajectory(filename, comment, self.best_path)

    def writeToTrajectory(self, filename, comment):
        self._writeToTrajectory(filename, comment, self.path)

    def _writeToTrajectory(self, filename, comment, path):
        trajectory = Trajectory(self.universe, filename, "w", comment)
        snapshot = SnapshotGenerator(self.universe,
                                     actions = [TrajectoryOutput(trajectory,
                                                                 ["all"],
                                                                 0, None, 1)])
        for step in path:
            self.universe.setConfiguration(step.conf)
            snapshot()
        trajectory.close()

    def initialize(self):
        if len(self.path[0]) > len(self.path[1]):
            move_from = 0
            move_to = 1
        else:
            move_from = 0
            move_to = 1
        while 1:
            distance = \
                 self.path[move_from][-1].distance(self.path[move_to][-1])
            if distance <= self.step_length: break
            start = self.path[move_from][-1]
            end = self.path[move_to][-1]
            direction = start.modeWeighting(start.differenceVectorTo(end))
            direction = start.differenceVectorTo(end)
            new = start.moveBy(direction, self.step_length)
            self.path[move_from].append(new)
            move_from, move_to = move_to, move_from
        self.path[1].reverse()
        self.path = self.path[0] + self.path[1]

    def lengthAndEquidistanceTerms(self):
        n = len(self.path)-1
        l = []
        for i in range(n):
            l.append(self.path[i].distance(self.path[i+1]))
        ln = max(0., N.sum(l)-n*self.step_length)
        ed = N.sum((N.array(l)-N.sum(l)/n)**2)
        return ln, ed

    def initializePenalties(self):
        n = len(self.path)-1
        self.global_term = N.zeros((n,), N.Float)
        self.deformation_term = N.zeros((n,), N.Float)
        for i in range(n):
            self.global_term[i], self.deformation_term[i] = \
                         self.path[i].stepPenalties(self.path[i+1])

    def updatePenalties(self, changed_step):
        i = changed_step
        self.global_term[i-1], self.deformation_term[i-1] = \
                             self.path[i-1].stepPenalties(self.path[i])
        self.global_term[i], self.deformation_term[i] = \
                             self.path[i].stepPenalties(self.path[i+1])
        
    def penalty(self):
        ln, ed = self.lengthAndEquidistanceTerms()
        ed = 2.e3*ed
        gl = 2.e4*N.sum(self.global_term)
        de = N.sum(self.deformation_term)
        sum = ln + ed + gl + de
        return sum

    def refine(self, step_count, protocol = None):
        if len(self.path) == 2:
            self.initialize()
        self.initializePenalties()
        p = self.penalty()
        if not hasattr(self, 'best_path'):
            self.best_path = copy(self.path)
            self.best_penalty = p
        if len(self.path) == 2:
            return
        if protocol is not None:
            protocol = open(protocol, 'w')
        temperature = 3*p
        while 1:
            to_be_moved = randrange(1, len(self.path)-1)
            displacement = randomParticleVector(self.universe, 1.)
            displacement = self.path[to_be_moved].modeWeighting(displacement)
            magnitude = gaussian(0., 0.1*self.step_length)
            old = self.path[to_be_moved]
            new = old.moveBy(displacement, magnitude)
            self.path[to_be_moved] = new
            self.updatePenalties(to_be_moved)
            new_p = self.penalty()
            if new_p > p and N.exp((p-new_p)/temperature) > uniform(0., 1.):
                self.path[to_be_moved] = old
                self.updatePenalties(to_be_moved)
            else:
                p = new_p
                if p < self.best_penalty:
                    self.best_penalty = p
                    self.best_path = copy(self.path)
            if protocol is not None:
                protocol.write('%f\n' % p)
                protocol.flush()
            step_count = step_count - 1
            if step_count == 0:
                break
        if protocol is not None:
            protocol.close()

def run(pdb1, pdb2, trajectory, nsteps, delpdb=0):
    universe = TransitionPathUniverse(CalphaForceField(2.5))
    universe.protein = Protein(pdb1, model='calpha')
    conf1 = copy(universe.configuration())
    struct2 = PDBConfiguration(pdb2).peptide_chains
    for i in range(len(universe.protein)):
        struct2[i].applyTo(universe.protein[i])
    if delpdb:
        os.unlink(pdb1)
        os.unlink(pdb2)
    tr, rms = universe.findTransformation(conf1)
    universe.applyTransformation(tr)
    conf2 =  copy(universe.configuration())
    universe.setBoundingBox(conf1, conf2)
    path = TransitionPath(universe, conf1, conf2, step_length=0.05, nmodes=50)
    path.refine(nsteps)
    path.writeBestToTrajectory(trajectory,
                               ("Transition path after %d steps, " % nsteps) +
                               ("energy %d," % path.best_penalty))
