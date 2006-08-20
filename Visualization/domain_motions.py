# Visualization of domain motions from DomainFinder output.
#
# Section and equation numbers in the comments refer to the article
#
# K. Hinsen, A. Thomas, M.J. Field
# Analysis of domain motions in large proteins,
# Proteins 34(3), 369-382 (1999)
#
from MMTK import *
from MMTK.Geometry import Line
from MMTK.Random import randomDirection
from Scientific.Geometry import Tensor, isVector
import Numeric, LinearAlgebra

# Load the data from DomainFinder. It is a dictionary containing the
# following entries:
#
# universe                  the universe containing the protein
# domains                   a list of domains (see below)
# conformation              the name of the PDB input file
# deformation               the deformation energies (a ParticleScalar object)
# deformation_threshold     the deformation threshold (a number)
# domain_coarseness         the domain coarseness (a number)
#
# If the domains were obtained from normal modes, there are two additional
# entries:
#
# nmodes                    the number of normal modes that were calculated
# mode_selection            the modes that were selected for domain analysis
#
# If the domains were obtained by comparing two conformations, there are
# also two additional entries:
#
# filter                    the noise filter setting (a number)
# comparison                the name of the second PDB input file
#
data = load('lyso.domains') # load output from DomainFinder
mode_number = 0             # index of the mode whose motions are shown
show_protein = 1            # Include a C-alpha model of the protein with
                            # domains shown in color. Set to zero if you
                            # plan to use a PDB file.
scale_factor = 500.         # to enhance visibility of domain motion arrows
similarity_threshold = 35.  # minimal "quality" of domains, empirical
pairs = [(1, 0)]            # domain pairs for which to show relative motions

# The domain list contains one Collection object for each domain.
# These collection objects each have two special attributes that describe
# the domains:
#
# similarity    an array containing similarity measures (Eq. 11)
#               for pairs of cubic subregions in the domain, or
#               the single number 0. if there is only one subregion.
#               The first value is the smallest similarity in the domain,
#               the others are the values for all other subregions
#               relative to one of the two most similar subregions,
#               sorted in ascending order. Therefore the number of
#               values is one smaller than the number of subregions.
#               This information can be used to decide how well-defined
#               a domain is.
# rbdata        an array containing the rigid-body motion parameters
#               for the cubic subregions. The first dimension of the
#               array enumerates the subregions. The second dimension
#               is of length 6*N, where N is the number of modes used
#               in the domain analysis (N=1 for configuration comparison).
#               The six parameters are the values of T and Phi as they
#               appear in Eqs. 10 and 11.
#
domains = data['domains']
universe = data['universe']
finite_displacements = not data.has_key('nmodes')
if finite_displacements:
    mode_number = 0

# Define the colors. First we set everything to white, then the
# atoms that belong to a domain get assigned a domain-specific color.
# We take into account only "good" domains, i.e. those whose
# similarity measure is below the user-defined similarity threshold.
colors = ['red', 'green', 'yellow', 'blue', 'orange', 'brown', 'violet',
          'magenta', 'cyan', 'olive']
for atom in universe.atomList():
    atom.color = 'white'
color_index = 0
for domain in domains:
    try:
        well_defined = domain.similarity[0] < similarity_threshold
    except TypeError:
        well_defined = 0
    if well_defined:
        domain.color = colors[color_index]
        for atom in domain.atomList():
            atom.color = colors[color_index]
        color_index = (color_index + 1) % len(colors)

# Import the graphics module. Here we use VMD, but you can use
# VRML or VRML2 as well. Just change these two lines to do so.
from Scientific.Visualization import VMD
graphics = VMD

# Define how to draw a spiral.
class Spiral(graphics.Group):
    
    """Spiral

    Constructor: Spiral(|point1|, |point2|, |start|,
                        |ratio|, |point_distance|,
                        |cylinder_radius|=0., **|attributes|)

    Arguments:

    |point1|, |point2| -- the end points of the spiral axis (vectors)

    |start| -- the first point of the spiral (a vector) or the radius
               (a positive number), which implies a randomly chosen
               start point

    |ratio| -- the ratio of movement along the axis to movement perpendicular
               to the axis (0 = circle, 1 = straight line)

    |point_distance| -- the distance between two neighbouring
                        points along the spiral

    |cylinder_radius| -- the radius of the cylinders of which the
                         spiral is made up. The default value is 0.,
                         which means that lines are used instead of spirals.

    |attributes| -- any graphics object attribute
    """
    
    def __init__(self, point1, point2, start, ratio,
                 point_distance, cylinder_radius=0., **attr):
        axis = (point2-point1).normal()
        angle = Numeric.arctan(ratio)
        if isVector(start):
            rvect = start-point1
            rvect = rvect-(rvect*axis)*axis
            radius = rvect.length()
        else:
            radius = start
            while 1:
                direction = axis.cross(randomDirection())
                if direction.length() > 1.e-3:
                    break
            start = point1 + radius*direction.normal()
        x = point_distance*Numeric.cos(angle)
        y = point_distance*Numeric.sin(angle)
        npoints = (point2-point1).length()/y
        tr = Translation(y*axis+point1) \
             * Rotation(axis, x/(2.*Numeric.pi*radius)) * Translation(-point1)
        p = start
        self.start_point = p
        objects = []
        for i in range(npoints):
            np = tr(p)
            if cylinder_radius == 0.:
                objects.append(apply(graphics.Line, (p, np), attr))
            else:
                objects.append(apply(graphics.Cylinder,
                                     (p, np, cylinder_radius), attr))
            p = np
        self.end_point = p
        graphics.Group.__init__(self, objects)



# The graphical representation of the protein consists of:
# 1) A colored sphere for each atom.
# 2) Grey cylinders joining the spheres.
protein = universe[0]   # There is nothing else in the universe.
scene = graphics.Scene(scale=1./Units.Ang)
if show_protein:
    for chain in protein:
        for residue in chain:
            atom = residue.C_alpha
            material = graphics.DiffuseMaterial(atom.color)
            scene.addObject(graphics.Sphere(atom.position(), 0.05,
                                            material = material))
        material = graphics.DiffuseMaterial('grey')
        for i in range(len(chain)-1):
            atom1 = chain[i].C_alpha
            atom2 = chain[i+1].C_alpha
            scene.addObject(graphics.Cylinder(atom1.position(),
                                              atom2.position(),
                                              0.01,
                                              material = material))

# For each domain, we indicate its absolute rigid-body motion by
# a cylinder representing the axis of a screw-motion description,
# plus a spiral around it that describes the screw motion.
for di1, di2 in pairs:
    domain1 = domains[di1]
    domain2 = domains[di2]
    if not hasattr(domain1, 'color'):
        print "Domain ", di1, " is not well defined"
        continue
    if not hasattr(domain2, 'color'):
        print "Domain ", di2, " is not well defined"
        continue
    parameters1 = Numeric.sum(domain1.rbdata)/len(domain1.rbdata)
    t1 = Vector(parameters1[mode_number:mode_number+3])
    phi1 = Vector(parameters1[mode_number+3:mode_number+6])
    parameters2 = Numeric.sum(domain2.rbdata)/len(domain2.rbdata)
    t2 = Vector(parameters2[mode_number:mode_number+3])
    phi2 = Vector(parameters2[mode_number+3:mode_number+6])
    # Use Eqs. 13 and 14.
    if finite_displacements: # conformation comparison
        rotation1 = Rotation(phi1.normal(), phi1.length())
        total1 = Translation(t1)*rotation1
        rotation2 = Rotation(phi2.normal(), phi2.length())
        total2 = Translation(t2)*rotation2
        relative = total1*total2.inverse()
        r_ref, rotation_axis, rotation_angle, translation = \
               relative.screwMotion()
    else: # normal modes
        phi = phi2-phi1
        rotation_axis = phi.normal()
        rotation_angle = phi.length()
        t = t2-t1
        translation = t*rotation_axis
        r_ref = phi.cross(t)/(phi*phi)
    # Find the point on the axis which is closest to the center of mass.
    r = Line(r_ref, rotation_axis).projectionOf(domain.centerOfMass())
    # Enhance the motion by a scale factor to make it better visible.
    translation = translation*scale_factor
    rotation_angle = rotation_angle*scale_factor
    # Add a cylinder and a spiral.
    material = graphics.DiffuseMaterial(domain1.color)
    scene.addObject(graphics.Cylinder(r-1.*Units.nm*rotation_axis,
                                      r+1.*Units.nm*rotation_axis,
                                      0.05, material = material))
    spiral_radius = 0.2*Units.nm
    line_length = Numeric.sqrt((spiral_radius*rotation_angle)**2
                               + translation**2)
    if abs(translation) < 0.001*spiral_radius*abs(rotation_angle):
        ratio = 1.
        segment_length = line_length
    else:
        ratio = abs(spiral_radius*rotation_angle/translation)
        segment_length = min(line_length/5., 0.05*Units.nm)
    scene.addObject(Spiral(r - 0.5*translation*rotation_axis,
                           r + 0.5*translation*rotation_axis,
                           spiral_radius, ratio, segment_length,
                           0.03, material = material))

# Show the complete scene
scene.view()
