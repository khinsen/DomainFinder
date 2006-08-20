import MMTK
import Numeric; N = Numeric
import LinearAlgebra; LA = LinearAlgebra
from Scientific.Geometry import Vector
from Scientific.Geometry.TensorModule import delta, epsilon


def rigidMovement(atoms, vector):
    a = N.zeros((len(atoms), 3, 2, 3), N.Float)
    b = N.zeros((len(atoms), 3), N.Float)
    for i in range(len(atoms)):
	a[i,:,0,:] = delta.array
	a[i,:,1,:] = (epsilon*atoms[i].position()).array
	b[i] = vector[atoms[i]].array
    a.shape = (3*len(atoms), 6)
    b.shape = (3*len(atoms),)
    vo = N.dot(LA.generalized_inverse(a), b)
    return Vector(vo[:3]), Vector(vo[3:])


def rigidRegions(vectors, deformation, cutoff,
		 cutoff_increment=1., subset=None):
    world = vectors[0].universe
    natoms = world.numberOfCartesianCoordinates()
    if subset is None:
	subset = world

    pc = MMTK.PartitionedAtomCollection(1.2, subset)
    regions = filter(lambda p: len(p) > 2,
		     map(lambda p: p[2], pc.partitions()))

    data = []
    rigid_regions = []
    for region in regions:
	r_data = []
	in_domain = 0
	flimit = cutoff
	for v, f in map(None, vectors, deformation):
	    f_mean = N.add.reduce(map(lambda r, f=f: f[r], region))/len(region)
	    in_domain = in_domain + (f_mean < flimit)
	    u, o = rigidMovement(region, v)
	    r_data = r_data + list(u) + list(o)
	    flimit = flimit*cutoff_increment
	if in_domain >= 0.5*len(vectors):
	    data.append(r_data)
	    rigid_regions.append(MMTK.Collection(region))

    data = N.array(data)
    arr = N.array(len(rigid_regions)*[None])
    for i in range(len(arr)):
	arr[i] = rigid_regions[i]

    return data, arr


def rigidRegionsFinite(comp, df, cutoff, subset=None):
    world = comp.universe
    natoms = world.numberOfCartesianCoordinates()
    if subset is None:
	subset = world

    pc = MMTK.PartitionedAtomCollection(1.2, world)
    regions = filter(lambda p: len(p) > 2, map(lambda p: p[2], pc.partitions()))

    data = []
    rigid_regions = []
    for region in regions:
	in_domain = 1
	df_mean = N.add.reduce(map(lambda a, f=df: f[a], region))/len(region)
	if df_mean < cutoff:	
	    r = MMTK.Collection(region)
	    tr, rms = r.findTransformation(comp)
	    v = tr.translation().displacement()
	    axis, angle = tr.rotation().axisAndAngle()
	    data.append(list(v)+list(angle*axis))
	    rigid_regions.append(r)

    data = N.array(data)
    arr = N.array(len(rigid_regions)*[None])
    for i in range(len(arr)):
	arr[i] = rigid_regions[i]

    return data, arr


def clusterAnalysis(matrix, tolerance = 10.):
    n = matrix.shape[0]
    clusters = []
    covered = N.zeros((n,), N.Int)
    while 1:
	pick = N.maximum.reduce(matrix)
	if N.maximum.reduce(pick) == 0.:
	    break
	imax = N.argmax(pick)
        max_sim = N.maximum.reduce(matrix[imax])
	cutoff = max_sim/tolerance
	mask = N.greater(matrix[imax], cutoff)
	similarity = N.sort(N.repeat(matrix[imax], mask))
	mask[imax] = 1
	clusters.append((similarity, mask))
        covered = N.logical_or(covered, mask)
	mask = N.logical_not(mask)
	matrix = matrix*(mask[N.NewAxis,:]*mask[:,N.NewAxis])
    if N.add.reduce(covered) < n:
	clusters.append((0., N.logical_not(covered)))
    return clusters
