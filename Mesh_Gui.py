# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 15:16:58 2020

@author: bluejgw
"""
import pyvista as pv
import sympy as sp
from sympy import Matrix, lambdify
import numpy as np
from PyQt5 import Qt, QtWidgets
from PyQt5.QtWidgets import QMessageBox
from pyvistaqt import QtInteractor
import sys, os
#import meshio
import trimesh

# from CGAL import CGAL_Polygon_mesh_processing
# current conda cgal is version 5.0.1, it doesn't include centroid()
# either wait till 5.0.3 is released on conda or DIY

# initiate stored mesh
mesh = pv.PolyData()

class MainWindow(Qt.QMainWindow):
    def __init__(self, parent=None, show=True):
        Qt.QMainWindow.__init__(self, parent)

        # create the frame
        self.frame = Qt.QFrame()
        vlayout = Qt.QVBoxLayout()

        # add the pyvista interactor object
        self.plotter = QtInteractor(self.frame)
        vlayout.addWidget(self.plotter.interactor)

        self.frame.setLayout(vlayout)
        self.setCentralWidget(self.frame)

        # simple menu
        mainMenu = self.menuBar()
        fileMenu = mainMenu.addMenu('File')
        editMenu = mainMenu.addMenu('Edit')
        
        # opening a mesh file
        self.open_mesh_action = Qt.QAction('Open Mesh...', self)
        self.open_mesh_action.triggered.connect(self.open_mesh)
        fileMenu.addAction(self.open_mesh_action)
        
        # exit button
        exitButton = Qt.QAction('Exit', self)
        exitButton.setShortcut('Ctrl+Q')
        exitButton.triggered.connect(self.close)
        fileMenu.addAction(exitButton)

        # inserting maximally inscribed cube
        self.max_cube_action = Qt.QAction('Max Cube', self)
        #self.max_cube_action.setCheckable(True)
        self.max_cube_action.triggered.connect(self.max_cube)
        editMenu.addAction(self.max_cube_action)
        
        #Create Cone in Mesh
        self.next_cubes_action = Qt.QAction('Next Cubes', self)
        self.next_cubes_action.triggered.connect(self.next_cubes)
        editMenu.addAction(self.next_cubes_action)
        
        # slice mesh horizontally based on internal cubes
        self.cube_hslice_action = Qt.QAction('Cube H-Slice', self)
        self.cube_hslice_action.triggered.connect(self.cube_hslice)
        editMenu.addAction(self.cube_hslice_action)

        # slice mesh (interactively)
        self.slice_action = Qt.QAction('Slice', self)
        self.slice_action.triggered.connect(self.slice)
        editMenu.addAction(self.slice_action)
        
        # slice mesh with clipping (interactively)
        self.clip_slice_action = Qt.QAction('Clip Slice', self)
        self.clip_slice_action.triggered.connect(self.clip_slice)
        editMenu.addAction(self.clip_slice_action)

        # create bounding box(es) for mesh (interactively)
        self.bounding_action = Qt.QAction('Bounding', self)
        self.bounding_action.triggered.connect(self.bounding_bar)
        editMenu.addAction(self.bounding_action)
        
        if show:
            self.show()

        self.plotter.add_axes(interactive=None, line_width=2, color=None, x_color=None, y_color=None, z_color=None, xlabel='X', ylabel='Y', zlabel='Z', labels_off=False, box=None, box_args=None)

    def open_mesh(self):
        """ add a mesh to the pyqt frame """
        global mesh

        # open file
        file_info = QtWidgets.QFileDialog.getOpenFileName()
        file_dir = file_info[0]
        
        # determine file type and if conversion needed
        head, tail = os.path.split(file_dir)
        root, ext = os.path.splitext(tail)

        # convert mesh file type
        #if ext != ".vtk" or ext != ".VTK":
        #    mesh = meshio.read(file_dir)
        #    meshio.write(root + ".vtk", mesh)
        #    mesh = pv.read(head + "/" + root + ".vtk")
            # need to store elsewhere or delete .vtk file in the future
        #else:
        #    mesh = pv.read(file_dir)

        # read mesh & transform according to principal axes
        pre = trimesh.load_mesh(file_dir)
        orient = pre.principal_inertia_transform
        pre = pre.apply_transform(orient)
        pre.export('data/'+ root + '_oriented.STL')
        mesh = pv.read('data/'+ root + '_oriented.STL')

        # show transformed mesh
        self.plotter.add_mesh(mesh, show_edges=True, color="w", opacity=0.6)

        # reset plotter
        self.reset_plotter()

        # find mesh centroid and translate the mesh so that's the origin
        self.centroid()

        # show bounding box
        self.plotter.add_bounding_box(opacity=0.5, color="y")

    def reset_plotter(self):
        """ clear plotter of mesh or interactive options """
        # clear plotter
        self.plotter.clear()
        self.plotter.clear_plane_widgets()
        self.plotter.reset_camera()
        self.update()
        
        # callback opened mesh
        self.plotter.add_mesh(mesh, show_edges=True, color="w", opacity=0.6)
        
        # show origin
        #self.plotter.add_axes_at_origin(xlabel='X', ylabel='Y', zlabel='Z', line_width=1, labels_off=False)

        # show floors
        #self.plotter.add_floor('-y')
        #self.plotter.add_floor('-z')
        
    def centroid(self):
        """ find centroid volumetrically and indicate on graph """
        global Vol_centroid, V

        # find the vertices & the vertex indices of each triangular face
        V = np.array(mesh.points)
        col = len(V)
        f_ind = np.array(mesh.faces.reshape((-1,4))[:, 1:4])
        
        # define an arbitrary start point from middle of max and min of X,Y,Z of
        # all points: in a convex manifold it falls inside the volume (requires
        # segmentation for general application)
        start = np.array(mesh.center)
        X_start = start[0]
        Y_start = start[1]
        Z_start = start[2]
        
        # initialize variables
        centroids = []
        Vol_total = 0
        Sum_vol_x = 0
        Sum_vol_y = 0
        Sum_vol_z = 0
        
        # find centroid from all tetrahedra made with arbitrary center and triangular faces
        for i in range(0, col-1, 3):          
            # find the center of each tetrahedron (average of X,Y,Z of 
            # 4 vertices, 3 from the triangle, and one arbitrary start point)
            X_cent = (X_start + V[f_ind[i,0],0] + V[f_ind[i+1,0],0] + V[f_ind[i+2,0],0]) / 4
            Y_cent = (Y_start + V[f_ind[i,1],1] + V[f_ind[i+1,1],1] + V[f_ind[i+2,1],1]) / 4
            Z_cent = (Z_start + V[f_ind[i,2],2] + V[f_ind[i+1,2],2] + V[f_ind[i+2,2],2]) / 4
    
            # compute the volume of each tetrahedron
            V1 = np.array([V[f_ind[i,0],0], V[f_ind[i,1],1], V[f_ind[i,2],2]])**2 - np.array([X_start, Y_start, Z_start])**2
            V2 = np.array([V[f_ind[i+1,0],0], V[f_ind[i+1,1],1], V[f_ind[i+1,2],2]])**2 - np.array([V[f_ind[i,0],0], V[f_ind[i,1],1], V[f_ind[i,2],2]])**2
            V3 = np.array([V[f_ind[i+2,0],0], V[f_ind[i+2,1],1], V[f_ind[i+2,2],2]])**2 - np.array([V[f_ind[i+1,0],0], V[f_ind[i+1,1],1], V[f_ind[i+1,2],2]])**2
            V1 = V1.reshape((-1,1))
            V2 = V2.reshape((-1,1))
            V3 = V3.reshape((-1,1))
            Vol = abs(np.linalg.det(np.hstack([V1, V2, V3]))) / 6
    
            # tally up each cycle
            Vol_total = Vol_total + Vol
            Sum_vol_x = Sum_vol_x + Vol * X_cent
            Sum_vol_y = Sum_vol_y + Vol * Y_cent
            Sum_vol_z = Sum_vol_z + Vol * Z_cent
            centroids.append([X_cent,Y_cent,Z_cent])
        
        # find & show centroid
        centroids = np.asarray(centroids)
        Vol_centroid = [Sum_vol_x, Sum_vol_y, Sum_vol_z] / Vol_total
        print("Total Volume:", Vol_total)
        print("Centroid:", Vol_centroid)

    def max_cube(self):
        """ add a maximally inscribed cube within the opened mesh """
        global max_c1, max_c2

        # reset plotter
        self.reset_plotter()

        # project cones to from centroid to find maximally inscribed cube
        ranges = mesh.bounds
        h = (abs(ranges[4]) + abs(ranges[5]))/3
        ang = np.arctan(0.5/(np.sqrt(2)/2))
        ang = float(90 - np.degrees(ang))
        max_c1 = pv.Cone(center=Vol_centroid+[0,0,h/2], direction=[0.0, 0.0, -1.0], height=h, radius=None, capping=False, angle=ang, resolution=100)
        max_c2 = pv.Cone(center=Vol_centroid-[0,0,h/2], direction=[0.0, 0.0, 1.0], height=h, radius=None, capping=False, angle=ang, resolution=100)
        #self.plotter.add_mesh(max_c1,color="r", opacity=0.2)
        #self.plotter.add_mesh(max_c2,color="r", opacity=0.2)

        # show centroid
        self.plotter.add_mesh(pv.PolyData(Vol_centroid), color='r', point_size=20.0, render_points_as_spheres=True)
        
        # find the nearest possible max cube vertex
        top = self.nearest_pt(max_c1, Vol_centroid)
        bottom = self.nearest_pt(max_c2, Vol_centroid)
        if top[0] < bottom[0]:
            p = top[1]
            V = top[2]
        else:
            p = bottom[1]
            V = bottom[2]
        
        # create and show max cube
        self.create_cube(V[p,:], Vol_centroid)
        self.plotter.add_mesh(max_cube, show_edges=True, color="b", opacity=0.6)
        #self.plotter.add_mesh(cell_center, color="r", point_size=8.0, render_points_as_spheres=True)

        # re-assign V as points of mesh
        V = np.array(mesh.points)

    def nearest_pt(self, cone, starting_pt):
        """ find nearest vertex: for segmented convex manifold, a cube with volume centroid as 
        center and nearest vertex as cube vertex, it falls inside the volume """
        global vert, p, clip

        clip = mesh.clip_surface(cone, invert=True)
        #clip2 = mesh.clip_surface(cone, invert=False)
        #diff = clip1.boolean_difference(clip2)
        #vert = np.array(diff.points)
        #self.plotter.add_mesh(clip, opacity=0.6, show_edges=True, color="g")

        # find nearest point in the clipped mesh
        vert = np.array(clip.points)
        c = len(vert)
        dist = np.zeros(c)
        for i in range(0, c):
            dist[i] = np.sqrt((vert[i,0] - starting_pt[0])**2 + (vert[i,1] - starting_pt[1])**2
                            + (vert[i,2] - starting_pt[2])**2)
                
        # find index of the nearest point
        nearest = min(dist)
        p = np.where(dist == nearest)
        p = p[0].item()

        return nearest, p, vert
            
    def create_cube(self, vertex, centroid):
        ''' create cube from the nearest pt & centroid '''
        global face_center, max_cube   

        # find the other 7 vertices
        # 3 vertices can be found by rotating the first point 90 degrees 3 times around Z axis of centroid
        # 4 vertices can be found by translating the first four vertices twice the half edge
        # found from the distance times sin(pi/4)
        t = sp.Symbol('t')
        Rz_t = Matrix([[sp.cos(t), -sp.sin(t), 0],
            [sp.sin(t), sp.cos(t), 0],
            [0, 0, 1]])
        #Txy = Matrix([[1, 0, centroid[0]],
        #    [0, 1, centroid[1]],
        #    [0, 0, 1]])
        #T_x_y = Matrix([[1, 0, -centroid[0]],
        #    [0, 1, -centroid[1]],
        #    [0, 0, 1]])
        #Rz_Txy_t = Txy * Rz_t * T_x_y
        #Rz_Txy = lambdify(t, Rz_Txy_t)
        Rz = lambdify(t, Rz_t)
        
        # construct the array of the first 4 vertices
        V_a = np.array(vertex-[centroid[0], centroid[1], 0])
        a_2 = np.dot(Rz(np.pi/2), V_a.T).T + np.array([centroid[0], centroid[1], 0])
        a_3 = np.dot(Rz(np.pi), V_a.T).T + np.array([centroid[0], centroid[1], 0])
        a_4 = np.dot(Rz(3*np.pi/2), V_a.T).T + np.array([centroid[0], centroid[1], 0])
        V_a = np.array(vertex)
        cube_V_start = np.array([V_a, a_2, a_3, a_4])
        print(cube_V_start)
        
        # find the translation distance
        half_edge = np.ones((4,1)) * [[0, 0, np.sign(centroid[2]-vertex[2])]] * np.sqrt((vertex[0]-centroid[0])**2 + (vertex[1]-centroid[1])**2) * sp.sin(sp.pi/4)
        translate = np.asarray(2*half_edge, dtype=np.float64)

        # construct the cube
        cube_V_end = np.add(cube_V_start, translate)
        cube_V = np.vstack((cube_V_start, cube_V_end))

        # construct the 6 faces of the cube
        cube_F =np.hstack([[4,0,1,2,3],
                        [4,0,3,7,4],
                        [4,0,1,5,4],
                        [4,1,2,6,5],
                        [4,2,3,7,6],
                        [4,4,5,6,7]])

        max_cube = pv.PolyData(cube_V, cube_F)

        cell_center = max_cube.cell_centers()
        face_center = np.array(cell_center.points)
        #print("Center of Faces:",face_center)
        
    def next_cubes(self):
        ''' create cubes within the mesh from the face centers of the first cube'''
        hi = 12
        ang = np.arctan(1/(np.sqrt(2)/2))
        ang = float(90 - np.degrees(ang))
        create_cone = pv.Cone(center=face_center[0] + [0,0,hi/2], direction = [0.0,0.0,-1.0], height = hi, radius=None, resolution= 100, angle = ang, capping=False)
        create_cone2 = pv.Cone(center=face_center[1] + [0,0,hi/2], direction = [-1.0,-1.0,0.0], height = hi, radius=None, resolution= 100, angle = ang, capping=False)
        create_cone5 = pv.Cone(center=face_center[5] + [0,0,-hi/2], direction = [0,0,1], height = hi, radius=None, resolution= 100, angle = ang, capping=False)
        clipped = mesh.clip_surface(create_cone, invert=True)
        clipped5 = mesh.clip_surface(create_cone5, invert=True)
        cone_intersection5 = np.array(clipped5.points)
        #nearest_cone = min(cone_intersection5)
        
        print("Bottom Intersection:",cone_intersection5)
        #cone_center = create_cone.cell_centers()
        #cone_center_points = np.array(cone_center.points)
        self.plotter.add_mesh(create_cone, color="y", opacity=0.6)
        self.plotter.add_mesh(clipped, show_edges=True, color="r", opacity=0.6)
        self.plotter.add_mesh(create_cone2, show_edges=True, color="y", opacity=0.6)
        self.plotter.add_mesh(create_cone5, show_edges=True, color="y", opacity=0.6)
        self.plotter.add_mesh(clipped5, show_edges=True, color="r", opacity=0.6)
        #clip1 = mesh.clip_surface(create_cone, invert=True)
        #clip2 = mesh.clip_surface(max_c2, invert=True)
        #clip = [clip1, clip2]
        #self.plotter.add_mesh(clip1, opacity=0.6, show_edges=True, color="w")
        #self.plotter.add_mesh(clip[1], opacity=0.6, show_edges=True, color="g")
        #self.plotter.add_mesh(cone_center, color="r", point_size=8.0, render_points_as_spheres=True)
        #print("Cone Center:",cone_center_points)

    def cube_hslice(self):
        """ slice mesh horizontally based on internal cubes """
        # reset plotter
        self.reset_plotter()

        # create sliced parts
        part1 = mesh.clip_closed_surface('z', origin=face_center[0])
        part2_a = mesh.clip_closed_surface('-z', origin=face_center[0])
        part2 = part2_a.clip_closed_surface('z', origin=face_center[5])
        part3 = mesh.clip_closed_surface('-z', origin=face_center[5])

        # display sliced parts
        self.plotter.clear()
        self.plotter.add_mesh(max_cube, show_edges=True, color="b", opacity=0.6)
        self.plotter.add_mesh(pv.PolyData(Vol_centroid), color='r', point_size=20.0, render_points_as_spheres=True)
        self.plotter.add_mesh(part1, show_edges=True, color="r", opacity=0.4)
        self.plotter.add_mesh(part2, show_edges=True, color="w", opacity=0.4)
        self.plotter.add_mesh(part3, show_edges=True, color="g", opacity=0.4)

    def slice(self):
        """ slice the mesh interactively """
        # reset plotter
        self.reset_plotter()

        self.plotter.add_mesh_slice_orthogonal(mesh)
    
    def clip_slice(self):
        """ slice & clip the mesh interactively """     
        # reset plotter
        self.reset_plotter()

        self.plotter.add_mesh_clip_plane(mesh)

    def bounding(self, level):
        level = int(level)
        bound = mesh.obbTree
        bound.SetMaxLevel(10)
        bound.GenerateRepresentation(level, boxes)
        self.plotter.add_mesh(boxes, opacity=0.2, color="g")
        return

    def bounding_bar(self):
        """ show various levels of OBB (Oriented Bounding Box) interactively """  
        # initialize bounding boxes mesh
        global boxes
        boxes = pv.PolyData()

        # reset plotter
        self.reset_plotter()

        self.plotter.add_slider_widget(self.bounding, [0, 10], title='Level')

    def closeEvent(self, event):
        reply = QMessageBox.question(self, "Window Close", "Are you sure you want to quit program?",
                                     QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
        if reply == QMessageBox.Yes:
            event.accept()
        else:
            event.ignore()
        
if __name__ == '__main__':
    app = Qt.QApplication(sys.argv)
    window = MainWindow()
    window.show()
    window.setWindowTitle("Mesh Visualization")
    QtWidgets.QApplication.setQuitOnLastWindowClosed(True)
    sys.exit(app.exec_())
