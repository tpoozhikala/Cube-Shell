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
import sys
#import os, meshio

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
        
        # indicate centroid
        self.show_cent_action = Qt.QAction('Show Centroid', self)
        self.show_cent_action.triggered.connect(self.show_cent)
        editMenu.addAction(self.show_cent_action)

        # inserting maximally inscribed cube
        self.max_cube_action = Qt.QAction('Max Cube', self)
        self.max_cube_action.triggered.connect(self.max_cube)
        editMenu.addAction(self.max_cube_action)
        
        # indicate slice lines on mesh according to major planes
        self.ortho_action = Qt.QAction('Ortho', self)
        self.ortho_action.triggered.connect(self.ortho)
        editMenu.addAction(self.ortho_action)

        # slice mesh (interactively)
        self.slice_action = Qt.QAction('Slice', self)
        self.slice_action.triggered.connect(self.slice)
        editMenu.addAction(self.slice_action)
        
        # slice mesh with clipping (interactively)
        self.clip_slice_action = Qt.QAction('Clip Slice', self)
        self.clip_slice_action.triggered.connect(self.clip_slice)
        editMenu.addAction(self.clip_slice_action)
        
        if show:
            self.show()

    def open_mesh(self):
        """ add a mesh to the pyqt frame """
        global mesh

        # open file
        file_info = QtWidgets.QFileDialog.getOpenFileName()
        file_dir = file_info[0]
        
        # determine file type and if conversion needed
        # head, tail = os.path.split(file_dir)
        # root, ext = os.path.splitext(tail)
        #if ext != ".vtk" or ext != ".VTK":
        #    mesh = meshio.read(file_dir)
        #    meshio.write(root + ".vtk", mesh)
        #    mesh = pv.read(head + "/" + root + ".vtk")
            # need to store elsewhere or delete .vtk file in the future
        #else:
        #    mesh = pv.read(file_dir)

        # read mesh
        mesh = pv.read(file_dir)

        # reset plotter
        self.reset_plotter()

        # find mesh centroid
        self.centroid()
    
    def reset_plotter(self):
        """ clear plotter of mesh or interactive options """
        # clear plotter
        self.plotter.clear()
        self.plotter.reset_camera()
        
        # callback opened mesh
        self.plotter.add_mesh(mesh, show_edges=True, color="w", opacity=0.6)
        
        # show floors
        self.plotter.add_floor('-y')
        self.plotter.add_floor('-z')

    def centroid(self):
        """ find centroid volumetrically and indicate on graph """
        global Vol_centroid, V, col

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
            V1 = np.array([V[f_ind[i,0],0], V[f_ind[i,1],1], V[f_ind[i,2],2]]) - np.array([X_start, Y_start, Z_start])
            V2 = np.array([V[f_ind[i+1,0],0], V[f_ind[i+1,1],1], V[f_ind[i+1,2],2]]) - np.array([V[f_ind[i,0],0], V[f_ind[i,1],1], V[f_ind[i,2],2]])
            V3 = np.array([V[f_ind[i+2,0],0], V[f_ind[i+2,1],1], V[f_ind[i+2,2],2]]) - np.array([V[f_ind[i+1,0],0], V[f_ind[i+1,1],1], V[f_ind[i+1,2],2]])
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

    def show_cent(self):
        """ indicate centroid of the mesh """
        # reset plotter
        self.reset_plotter()

        self.plotter.add_mesh(pv.PolyData(Vol_centroid), color='r', point_size=20.0, render_points_as_spheres=True)

    def max_cube(self):
        """ add a maximally inscribed cube within the opened mesh """
        # reset plotter
        self.reset_plotter()

        # project cones to from centroid to find maximally inscribed cube
        h = 10
        ang = np.arctan(0.5/(np.sqrt(2)/2))
        ang = float(90 - np.degrees(ang))
        c1 = pv.Cone(center=Vol_centroid+[0,0,h/2], direction=[0.0, 0.0, -1.0], height=h, radius=None, capping=True, angle=ang, resolution=100)
        c2 = pv.Cone(center=Vol_centroid-[0,0,h/2], direction=[0.0, 0.0, 1.0], height=h, radius=None, capping=True, angle=ang, resolution=100)
        self.plotter.add_mesh(c1,color="g", opacity=0.6)
        self.plotter.add_mesh(c2,color="g", opacity=0.6)

        #inter = mesh.boolean_cut(c1)
        #self.plotter.add_mesh(inter,show_edges=False, color="w", opacity=0.6)

        # find nearest vertex: for segmented convex manifold, a cube with volume centroid as 
        # center and nearest vertex as cube vertex, it falls inside the volume
        dist = np.zeros(col-1)
        for i in range(0, col-1):
            dist[i] = np.sqrt((V[i,0] - Vol_centroid[0])**2 + (V[i,1] - Vol_centroid[1])**2
                              + (V[i,2]-Vol_centroid[2])**2)
        nearest = min(dist)
        
        # find index of the nearest vertex
        p = np.where(dist == nearest)
        p = p[0].item()
        
        # find the 7 other vertices
        # for axisymmetric parts
        # 3 vertices can be found by rotating the first point 90 degrees 3 times around Z axis
        # 4 vertices can be found by translating the first four twice the half edge
        # found from the distance times sin(pi/4)
        t = sp.Symbol('t')
        Rz_t = Matrix([[sp.cos(t), -sp.sin(t), 0],
              [sp.sin(t), sp.cos(t), 0],
              [0, 0, 1]])
        Rz = lambdify(t, Rz_t)
        
        V_a = np.array(V[p,:])
        a_2 = np.dot(Rz(np.pi/2), V_a.T).T
        a_3 = np.dot(Rz(np.pi), V_a.T).T
        a_4 = np.dot(Rz(3*np.pi/2), V_a.T).T
        cube_V_mid = np.array([V_a, a_2, a_3, a_4])
        
        half_edge = np.ones((4,1)) * [[0, 0, np.sign(Vol_centroid[2]-V[p,2])]] * np.sqrt(V[p,0]**2 + V[p,1]**2) * sp.sin(sp.pi/4)
        half_edge = np.asarray(half_edge, dtype=np.float64)

        cube_V_top = np.add(cube_V_mid, half_edge)
        cube_V_bottom = np.subtract(cube_V_mid, half_edge)
        cube_V = np.vstack((cube_V_top, cube_V_bottom))
        cube_F =np.hstack([[4,0,1,2,3],
                           [4,0,3,7,4],
                           [4,0,1,5,4],
                           [4,1,2,6,5],
                           [4,2,3,7,6],
                           [4,4,5,6,7]])

        max_cube = pv.PolyData(cube_V, cube_F)
        #self.plotter.add_mesh(max_cube, show_edges=True, color="b", opacity=0.6)

        cube_test =np.hstack([[4,0,1,2,3]])
        test = pv.PolyData(cube_V_mid,cube_test)
        self.plotter.add_mesh(test, show_edges=True, color="b", opacity=0.6)
        
    def ortho(self):
        """ indicate slice lines according to major planes """
        # reset plotter
        self.reset_plotter()

        slcOrtho = mesh.slice_orthogonal()
        self.plotter.add_mesh(slcOrtho, color="r")
    
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
