import bpy
import numpy as np
import time

stride_size_v = 8*3
stride_size_f = 6*4
def add_cube_vertex(index_0,pos,scales):
    v = []
    e = []
    v = np.array([[-1.0, -1.0, -1.0],
                  [-1.0, 1.0, -1.0],
                  [1.0, 1.0, -1.0],
                  [1.0, -1.0, -1.0],
                  [-1.0, -1.0, 1.0],
                  [-1.0, 1.0, 1.0],
                  [1.0, 1.0, 1.0],
                  [1.0, -1.0, 1.0]],dtype=np.float32)
    scale_x,scale_y,scale_z = scales
    pos_x,pos_y,pos_z = pos
    v[:,0] = (v[:,0] * scale_x) + pos_x
    v[:,1] = (v[:,1] * scale_y) + pos_y
    v[:,2] = (v[:,2] * scale_z) + pos_z
    faces = np.array([[4,5,1,0],
                      [5,6,2,1],
                      [6,7,3,2],
                      [7,4,0,3],
                      [0,1,2,3],
                      [7,6,5,4]],dtype=np.int32)+(8*i)
    return v,faces

def returnCubeMesh(passed_name):
    v = []
    e = []

    # Points that define the default cube.
    v = np.array([[-1.0, -1.0, -1.0],
                  [-1.0, 1.0, -1.0],
                  [1.0, 1.0, -1.0],
                  [1.0, -1.0, -1.0],
                  [-1.0, -1.0, 1.0],
                  [-1.0, 1.0, 1.0],
                  [1.0, 1.0, 1.0],
                  [1.0, -1.0, 1.0]],dtype=np.float32)
    # for vertex in vertices:
    #     v.append(Vector(vertex))
    # Default cube face order.
    faces = [[4,5,1,0],
             [5,6,2,1],
             [6,7,3,2],
             [7,4,0,3],
             [0,1,2,3],
             [7,6,5,4]]
    me = bpy.data.meshes.new(passed_name)
    me.from_pydata(v, e, faces)
    me.update()
    me.validate(verbose = True)  # useful for development when the mesh may be invalid.
    return me

def add_cube(pos,scales):
    scale_x,scale_y,scale_z = scales
    pos_x,pos_y,pos_z = pos
    mesh = returnCubeMesh("mpset")
    obj = bpy.data.objects.new("mpset",mesh)
    # obj.data = mesh
    return obj

objs = bpy.data.objects
scn = bpy.context.scene

# scene_objs = [
#         objs['mp']
#         ]
# original_mp = objs['mp']
csv_dir = bpy.path.abspath("//csvs\\")
voxel_data = np.loadtxt(csv_dir+"./frame_00090.csv",delimiter=",")
voxel_data = voxel_data[:,:]
my_coll = bpy.data.collections.new("MyCollection")
# bpy.ops.mesh.primitive_cube_add()
# ob = bpy.context.object
# copies = []
## add objects
start = time.time()
v = np.zeros((len(voxel_data)*8,3),dtype=np.float32)
e = np.array([],dtype=np.int32)
f = np.zeros((len(voxel_data)*6,4),dtype=np.int32)
for i,mp in enumerate(voxel_data):
    scale = 0.01
    mp[0:6] = mp[0:6]*scale
    mp[3:6] = mp[3:6]*0.5
    x,y,z = mp[0],mp[2],mp[1]
    lx,ly,lz = mp[3],mp[5],mp[4]
    # copy = ob.copy()
    # copy.data = ob.data.copy()
    # copy.location.x = i
    # bpy.context.collection.objects.link(copy)
    # obj = bpy.data.objects.new('mpset',None)
    # copy = add_cube((x,y,z),(lx,ly,lz))
    vi,fi = add_cube_vertex(i,(x,y,z),(lx,ly,lz))
    v[i*8:(i+1)*8,:] = vi
    f[i*6:(i+1)*6,:] = fi
    # v = np.concatenate((v,vi))
    # f = np.concatenate((f,fi))

me = bpy.data.meshes.new("meshdata")
me.from_pydata(v, [], f)
me.update()
#me.validate(verbose = True)  # useful for development when the mesh may be invalid.
me.validate()
obj = bpy.data.objects.new("mpset",me)
my_coll.objects.link(obj)
bpy.context.scene.collection.children.link(my_coll)

print(v)
print(f)
end = time.time()
print(end - start)


# degp = bpy.context.evaluated_depsgraph_get()
# object = bpy.data.objects["ps"]
# particle_systems = object.evaluated_get(degp).particle_systems
# particles = particle_systems[0].particles
# totalParticles = len(particles)
# flatList = [0]*(3*totalParticles)

# # additionally set the location of all particle locations to [0, 0, 0]
# particles.foreach_set("location", flatList)

# for collection in bpy.data.collections:
#    print(collection.name)
