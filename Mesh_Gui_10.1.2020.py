# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 15:16:58 2020

@author: bluej
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

#from CGAL import CGAL_Polygon_mesh_processing
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
        self.max_cube_action.triggered.connect(self.max_cube)
        editMenu.addAction(self.max_cube_action)
        
        # slice mesh
        self.slice_action = Qt.QAction('Slice', self)
        self.slice_action.triggered.connect(self.slice)
        editMenu.addAction(self.slice_action)
        
        if show:
            self.show()

    def open_mesh(self):
        """ add a mesh to the pyqt frame """
        self.plotter.clear()
        
        file_info = QtWidgets.QFileDialog.getOpenFileName()
        file_dir = file_info[0]
        
        global mesh
        
        # determine mesh type and if conversion needed
        head, tail = os.path.split(file_dir)
        root, ext = os.path.splitext(tail)
        #if ext != ".vtk" or ext != ".VTK":
        #    mesh = meshio.read(file_dir)
        #    meshio.write(root + ".vtk", mesh)
        #    mesh = pv.read(head + "/" + root + ".vtk")
            # need to store elsewhere or delete .vtk file in the future
        #else:
        #    mesh = pv.read(file_dir)
        mesh = pv.read(file_dir)
        
        #centers = mesh.cell_centers()
        
        self.plotter.add_mesh(mesh, show_edges=True, color="w", opacity=0.6)
        
        # show floors
        self.plotter.add_floor('-y')
        self.plotter.add_floor('-z')
        
        #self.plotter.add_mesh(mesh.extract_all_edges(), color="k", line_width=1)
        #self.plotter.add_mesh(centers, color="r", point_size=8.0, render_points_as_spheres=True)
        self.plotter.reset_camera()
    
    def max_cube(self):
        """ add a maximally inscribed cube within the opened mesh """
        V = np.array(mesh.points)
        col = len(V)
        #print(V)
        #print(col)
        
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
        
        # find centroid from all tetrahedra made with arbitrary center and STL
        # triangles
        for i in range(0, col-1, 3):
            # save the coordinates of the verteices of each triangle
            #X=[[X_start, X_start, X_start, V(i,0)],
            #   [V(i,0), V(i+1,0), V(i+2,0), V(i+1,0)],
            #   [V(i+1,0), V(i+2,0), V(i,0), V(i+2,0)]]
            #Y=[[Y_start, Y_start, Y_start, V(i,1)],
            #   [V(i,1), V(i+1,1), V(i+2,1), V(i+1,1)],
            #   [V(i+1,1), V(i+2,1), V(i,1), V(i+2,1)]]
            #Z=[[Z_start, Z_start, Z_start, V(i,2)],
            #   [V(i,2), V(i+1,2), V(i+2,2), V(i+1,2)],
            #   [V(i+1,2), V(i+2,2), V(i,2), V(i+2,2)]]
            
            # find the center of each tetrahedron (average of X,Y,Z of 
            # 4 vertices, 3 from the triangle, and one arbitrary start point)
            X_cent = (X_start + V[i,0] + V[i+1,0] + V[i+2,0]) / 4
            Y_cent = (Y_start + V[i,1] + V[i+1,1] + V[i+2,1]) / 4
            Z_cent = (Z_start + V[i,2] + V[i+1,2] + V[i+2,2]) / 4
    
            # compute the volume of each tetrahedron
            V1 = np.array([V[i,0], V[i,1], V[i,2]]) - np.array([X_start, Y_start, Z_start])
            V2 = np.array([V[i+1,0], V[i+1,1], V[i+1,2]]) - np.array([V[i,0], V[i,1], V[i,2]])
            V3 = np.array([V[i+2,0], V[i+2,1], V[i+2,2]]) - np.array([V[i+1,0], V[i+1,1], V[i+1,2]])
            V1 = V1.reshape((-1,1))
            V2 = V2.reshape((-1,1))
            V3 = V3.reshape((-1,1))
    
            Vol = abs(np.linalg.det(np.hstack([V1, V2, V3]))) / 6
            #print(Vol)
    
            Vol_total = Vol_total + Vol
            Sum_vol_x = Sum_vol_x + Vol * X_cent
            Sum_vol_y = Sum_vol_y + Vol * Y_cent
            Sum_vol_z = Sum_vol_z + Vol * Z_cent
    
            centroids.append([X_cent,Y_cent,Z_cent])
        centroids = np.asarray(centroids)
        
        #print(centroids)
        Vol_centroid = [Sum_vol_x, Sum_vol_y, Sum_vol_z] / Vol_total
        #print(Sum_vol_x)
        #print(Sum_vol_y)
        #print(Sum_vol_z)
        print("Total Volume:", Vol_total)
        print("Centroid:", Vol_centroid)
        
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
        #cube_V_top = np.array([V_a, a_2, a_3, a_4])
        cube_V_mid = np.array([V_a, a_2, a_3, a_4])
        
        half_edge = np.ones((4,1)) * [[0, 0, np.sign(Vol_centroid[2]-V[p,2])]] * np.sqrt(V[p,0]**2 + V[p,1]**2) * sp.sin(sp.pi/4)
        half_edge = np.asarray(half_edge, dtype=np.float64)
        #top_to_bottom = np.ones((4,1)) * [[0, 0, 2*np.sign(Vol_centroid[2]-V[p,2])]] * np.sqrt(V[p,0]**2 + V[p,1]**2) * sp.sin(sp.pi/4)
        #top_to_bottom = np.asarray(top_to_bottom, dtype=np.float64)

        cube_V_top = np.add(cube_V_mid, half_edge)
        cube_V_bottom = np.subtract(cube_V_mid, half_edge)
        #cube_V_bottom = np.add(cube_V_top, down)
        cube_V = np.vstack((cube_V_top, cube_V_bottom))
        
        cube_F =np.hstack([[4,0,1,2,3],
                           [4,0,3,7,4],
                           [4,0,1,5,4],
                           [4,1,2,6,5],
                           [4,2,3,7,6],
                           [4,4,5,6,7]])
        max_cube = pv.PolyData(cube_V, cube_F)
        self.plotter.add_mesh(max_cube, show_edges=True, color="b", opacity=0.6)
        #f_top = np.array([4,0,1,2,3])
        #top = pv.PolyData(cube_V_top, f_top)
        #self.plotter.add_mesh(top, show_edges=True, color="b", opacity=0.6)
        
    def slice(self):
        """ slice the mesh according to user input """
        slcOrtho = mesh.slice_orthogonal()
        self.plotter.add_mesh(slcOrtho, color="r")
        self.plotter.reset_camera()
    
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
