# coding=utf-8
# ------------------------------Documentation------------------------------------------------------
# This code generates different 2D material configuration with variate file format, such as LAMMPS.
# Current 2D material type: 1. Graphene 2. MoS2 3. Diamond Surface/nanopillar 4. Li Surface/nanopillar
# Current configuration file format: 1. LAMMPS
#                                    2. VESTA .xyz
# 
# Units: [L] = Angstrom
#        [E] = eV
#        [Angle] = deg
#        
# Last modified: 2/21/2019 by Changhao Li
# -------------------------------------------------------------------------------------------------
import os, sys
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import math
import copy

A0 = 2.02569 # Lattice Parameter for MoS2
Dx = 1.01284; Dy = 1.59481; Dz = 1.56534;

def search_word(FILE, sword, FIMP):
    ifs = open(FILE, "r")
    nline = sum(1 for line in open(FILE))
    for i in range(nline):
        line = ifs.readline()
        data = line.split()
        if len(data) == 0:
            continue
        if line[0] == "#":
            continue
        if data[0] == sword:
            if data[0] == "Norbital":
                ndat = len(data) - 1
                Norb = np.zeros(ndat, dtype=int)
                for i in range(ndat):
                    Norb[i] = int(data[1+i])
                return Norb
            else:
                ifs.close()
                return data[1]
    if FIMP == 1:
        print("Error: Cannot find %s in file %s"%(sword, FILE))
        sys.exit()
    return 0

def file_check(fname):
    if os.path.exists(fname) == False:
        print("Error: %s does not exist"%fname)
        sys.exit()
    return 0
 
def IsInside(x1,y1,x2,y2,x3,y3,x,y):
    # if (x, y) is in the triangle composed of (x1, y1), (x2, y2), (x3, y3), return True
    def IsTrangleOrArea(x1,y1,x2,y2,x3,y3):
        return abs((x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2)) / 2.0)

    ABC = IsTrangleOrArea(x1,y1,x2,y2,x3,y3)
    PBC = IsTrangleOrArea(x,y,x2,y2,x3,y3)
    PAC = IsTrangleOrArea(x1,y1,x,y,x3,y3)
    PAB = IsTrangleOrArea(x1,y1,x2,y2,x,y)
 
    return (ABC == PBC + PAC + PAB)

def deg2rad(ang):
    return ang/180.0 * math.pi

class InputParams:
    def __init__(self, fname):
        self.matType = ""
        self.phase = ""
        self.coordMethod = ""
        self.outputStyle = ""
        self.read_input(fname)

    def read_input(self, fname):
        file_check(fname)
        self.matType = search_word(fname, "matType", 0)
        self.phase = search_word(fname, "phase", 0)
        self.coordMethod = search_word(fname, "coordMethod", 0)
        self.outputStyle = search_word(fname, "outputStyle", 0)


class Atom:
    def __init__(self):
        self.x = ""
        self.y = ""
        self.z = ""
        self.t1 = ""
        self.t2 = ""
        self.t3 = ""
        self.type = 1
        self.num = 0
        self.DEL = 0

    def __init__(self, xx, yy, zz, TYPE, NUM):
        self.x = xx
        self.y = yy
        self.z = zz
        self.type = TYPE
        self.num = NUM

    def translate(self, xx, yy, zz):
        self.x += xx
        self.y += yy
        self.z += zz

    def printCord(self):
        print("Atom %s: %s, x = %s, y = %s, z = %s\n"%(self.num, self.type, self.x, self.y, self.z))

class Material2D:
    def __init__(self):
        self.natom = 0
        self.atomList = [];
        self.ntype = 0
        self.box = []
        self.coordMethod = "userDefined"
        self.outputStyle = "LAMMPS"

    def readFromXYZ(self, fname):
        a = np.array([0.0, 0.0, 0.0]); b = a.copy(); c = a.copy()
        nline = 0
        with open(fname, 'r') as f:
            content = f.readlines()
            for line in content:
                nline += 1
                line = line.strip().split()
                if len(line) >= 2:
                    if line[0] == 'Position':
                        break
                    if line[-1] == 'atoms':
                        self.natom = int(line[0])
                    if line[-2] == 'xlo' and line[-1] == 'xhi':
                        a[0] = float(line[1]) - float(line[0])
                    if line[-2] == 'ylo' and line[-1] == 'yhi':
                        b[1] = float(line[1]) - float(line[0])
                    if line[-2] == 'zlo' and line[-1] == 'zhi':
                        c[2] = float(line[1]) - float(line[0])
                    if len(line) >= 3:
                        if line[-3] == 'xy':
                            b[0] = float(line[0]); c[0] = float(line[1]); c[1] = float(line[2])
            self.box = np.vstack((a, b, c))

            for i in range(self.natom):
                line = content[nline + i]
                line = line.strip().split()
                iatom = Atom(float(line[3]), float(line[4]), float(line[5]), int(line[1]), int(line[0]))
                self.atomList.append(iatom)

    def triclincLattice(self, l, m, n):
        unit = Material2D()
        unit.copy(self)
        p = 0
        for i in range(l):
            for j in range(m):
                for k in range(n):
                    if i == 0 and j == 0 and k == 0:
                        continue                  
                    t = i*self.box[0, :] + j*self.box[1, :] + k*self.box[2, :]
                    iunit = unit.translate_copy(t[0], t[1], t[2])
                    self.merge(iunit)
                    p += 1
        self.box = np.vstack((l*self.box[0, :], m*self.box[1, :], n*self.box[2, :]))
        print("%s lattices in total"%(p + 1))

    def output_triclinc(self, ofilename, structName, charge=False, isTriclinc=True):
        with open(ofilename, 'w') as f:
            f.write("Configuration for %s Bulk Crystal\n"%structName)
            f.write("\n")
            f.write("%s atoms\n"%len(self.atomList))
            f.write("1 atom types\n")
            if isTriclinc:
                f.write("ITEM: BOX BOUNDS xy xz yz xx yy zz\n")
            f.write("%s %s xlo xhi\n"%(0, self.box[0, 0]))
            f.write("%s %s ylo yhi\n"%(0, self.box[1, 1]))
            f.write("%s %s zlo zhi\n"%(0, self.box[2, 2]))
            if isTriclinc:
                f.write("%s %s %s xy xz yz\n"%(self.box[1, 0], self.box[2, 0], self.box[2, 1]))
            f.write("\n")
            f.write("Atoms\n")
            f.write("Position for %s atoms\n"%structName)
            if charge == False:
                for i in range(len(self.atomList)):
                    iatom = self.atomList[i]
                    f.write("%s %s %12.6f %12.6f %12.6f\n"%(i+1, iatom.type, iatom.x, iatom.y, iatom.z))
            else:
                for i in range(len(self.atomList)):
                    iatom = self.atomList[i]
                    f.write("%s %s 0 %12.6f %12.6f %12.6f\n"%(i+1, iatom.type, iatom.x, iatom.y, iatom.z))


    def cut3D(self, rule, norm):
        new_list = []
        for atom in self.atomList:
            if rule(atom, norm):
                new_list.append(atom)
        self.atomList = new_list
        return self

    def cylinder_cut(self, xx, yy, r, zlo, zhi, atom, norm):
        if ((xx - atom.x)**2 + (yy - atom.y)**2 <= r**2) and (atom.z >= zlo) and (atom.z < zhi):
            if norm == 1:
                return 0
            elif norm == 0:
                return 1
            else:
                print("cylinder_cut: false norm value\n")

    def semi_sphere_cut(self, xx, yy, zz, r, atom, norm):
        if (xx - atom.x)**2 + (yy - atom.y)**2 + (zz - atom.z)**2 <= r**2 and atom.z >= zz:
            if norm == 1:
                return 0
            elif norm == 0:
                return 1
            else:
                print("semi_sphere_cut: false norm value\n")

    def surface_cut(self, xlo, xhi, ylo, yhi, zlo, zhi, atom, norm):
        if (atom.z >= zlo and atom.z < zhi) and \
           (atom.x >= xlo and atom.x < xhi) and \
           (atom.y >= ylo and atom.y < yhi):
            if norm == 1:
                return 0 
            elif norm == 0:
                return 1
            else:
                print("surface_cut: false norm value\n")


    def pure_bending_mapping(self, radius, length):
        r = radius; l = length
        for atom in self.atomList:
            x = atom.x; y = atom.y; z = atom.z
            theta = z / radius
            atom_r = r - atom.x
            atom.x = x + atom_r * (1 - math.cos(theta))
            atom.z = atom_r * math.sin(theta)
            


    def copy(self, mat2D):
        self.natom = mat2D.natom
        self.atomList = copy.deepcopy(mat2D.atomList)
        self.ntype = mat2D.ntype
        self.box = copy.deepcopy(mat2D.box)
        self.coordMethod = mat2D.coordMethod
        self.outputStyle = mat2D.outputStyle
        return self

    def print_info(self):
        print("natom: %s\natomList: %s\nntype: %s\n"%(self.natom, self.atomList, self.ntype))

    def make_configuration(self, myinput):
        if self.coordMethod == "userDefined":
            self.userDefined()

        if self.coordMethod == "rectangle":
            pass

    def hBN_unit_cell(self, C=1):
        # Return a list which contains a unit cell of hBN (2 Mo atoms, 4 S atoms, 2H phase)
        # Atom type: 1 -> S
        #            2 -> Mo
        # Dx, Dy, Dz, A0 have been defined in the beginning as global variables
        dx = Dx * C; dy = Dy * C; dz = Dz * C; A = A0 * C
        Mo1 = Atom(0, 0, 0, 2, 1) #Mo
        Mo2 = Atom(dx+A, dy, 0, 2, 2) #Mo
        S1 = Atom(dx, dy, dz, 1, 3) #S
        S2 = Atom(dx, dy, -dz, 1, 4) #S
        S3 = Atom(2*dx+A, 0, dz, 1, 5) #S
        S4 = Atom(2*dx+A, 0, -dz, 1, 6) #S
        return [Mo1, Mo2, S1, S2, S3, S4]

    def merge(self, mat2D):
        self.natom += mat2D.natom
        for atom in mat2D.atomList:
            self.atomList.append(atom)
        self.natom = len(self.atomList)
        return self

    def translate(self, xx, yy, zz=0):
        for atom in self.atomList:
            atom.x += xx
            atom.y += yy
            atom.z += zz
        # Here returns self, which has been translated
        return self

    def translate_copy(self, xx, yy, zz=0):
        cp = Material2D()
        cp.copy(self)
        for atom in cp.atomList:
            atom.x += xx
            atom.y += yy
            atom.z += zz
        # Here returns the copy of translated sheets, not self
        return cp

    def rotate(self, tz, ty=0, tx=0, cx=0, cy=0, cz=0):
        # Only Z-direction rotation is implemented
        # Unit: Degree
        tz = tz / 180.0 * math.pi; ty = ty / 180.0 * math.pi; tx = tx / 180.0 * math.pi
        rMatrix = np.array([[math.cos(tz), -math.sin(tz), 0],\
                            [math.sin(tz),  math.cos(tz), 0],\
                            [           0,             0, 1]])
        for atom in self.atomList:
            posVec = np.array([atom.x, atom.y, atom.z])
            centVec = np.array([cx, cy, cz])
            r = posVec - centVec
            posVec = np.dot(rMatrix, r.T) + centVec
            atom.x = posVec[0]; atom.y = posVec[1]; atom.z = posVec[2]

        return self

    def cut(self, x1, y1, x2, y2, norm):
        # Cut the sheet with the line whose vertices are (x1, y1) and (x2, y2). 
        # norm is the direction of cutting
        # cut_tol controls the edge
        new_list = [];
        cut_tol = 0;
        # The common cases, in which k is finite value
        if abs(x1 - x2) > 1e-6:
            k = (y1 - y2) / (x1 - x2)
            for atom in self.atomList:
                x = atom.x; y = atom.y
                if norm == 0:
                    if  (k * (x - x1)) >= (y - y1) or abs(k*(x-x1)-(y-y1)) < cut_tol:
                        new_list.append(atom)
                if norm == 1:
                    if (k * (x - x1)) < (y - y1) or abs(k*(x-x1)-(y-y1)) < cut_tol:
                        new_list.append(atom)
        # Degenerated case, in which k is infinite
        else:
            for atom in self.atomList:
                x = atom.x; y = atom.y
                if norm == 0:
                    if x >= x1:
                        new_list.append(atom)
                if norm == 1:
                    if x < x1:
                        new_list.append(atom)

        self.atomList = new_list

        return self

    def remove_redundant(self, cut_r):
        # Remove the atom whose coordinates are the same or too close to other atoms
        # cur_r controls the tolerance
        def isAtomInList(aList, a, cut_r):
            for atom in aList:
                if (a.x - atom.x)**2 + (a.y - atom.y)**2 + (a.z - atom.z)**2 <= cut_r**2:
                    return True
            return False
        def isAtomInSurface(Mat2D, atom, cut_r):
            if abs(atom.x - Mat2D.box[0, 0]) < cut_r or abs(atom.y - Mat2D.box[1, 1]) < cut_r or abs(atom.z - Mat2D.box[2, 2]) < cut_r:
                return True
            return False

        new_atom_list = []
        for iatom in self.atomList:
            if not (isAtomInList(new_atom_list, iatom, cut_r) or isAtomInSurface(self, iatom, cut_r)):
                new_atom_list.append(iatom)

        self.atomList = new_atom_list
        self.natom = len(new_atom_list)


    def rectangle_sheet(self, Tx, Ty, C=1):
        # Make a rectangle sheet of 2D materials with the repetation of (Tx, Ty) unit cells
        # The unit cell is defined in the method of Material2D.
        unit = Material2D()
        unit.natom = 6
        unit.atomList = unit.hBN_unit_cell()
        A = A0 * C; xx = 3 * A; yy = A * math.sqrt(3)
        for i in range(0, Tx):
            for j in range(0, Ty):
                self.merge(unit.translate_copy(xx*i, yy*j))

    def find_atom_in_cubic(self, xlo, xhi, ylo, yhi, zlo, zhi):
        # return a list which contains all the atom in the cubic 
        # The size of the cubic is defined by (xlo, xhi, ylo, yhi, zlo, zhi)
        pass

    def find_n_m(self, theta, tol=1e-3):
        # Given rotational angle theta, find the lattice vector (n, m)
        # The default tolerance is 1e-2
        # Unit: Angle
        theta = theta % 120.0
        theta = theta / 180.0 * math.pi
        a1 = 1 - math.sqrt(3)*math.tan(theta)
        b1 = -1 - math.sqrt(3)*math.tan(theta)
        a2 = 1
        b2 = 0
        coeff1 = np.array([[a1, b1], [a2, b2]])
        coeff2 = np.array([0, 1])
        x = np.linalg.solve(coeff1, coeff2)
        if False:
            n = x[0]; m = x[1]; i = 2
            while not(abs(m%1 - 1) <= tol or abs(m%1) <= tol):
                n *= i/(i-1); m *= i/(i-1); i += 1
        return x

# -------------------------- User-defined method --------------------------------------- #   
    def userDefined(self, C=1, method=1):
        # User-defined method to generate complex multi-phase 2D materials sheets
        # C controls the size (1 is default)
        print("generate configuration: user-defined method\n")
        A = A0 * C; xx = 3 * A; yy = A * math.sqrt(3)

        method = 1

        def sheet_cut_merge(sheet1, sheet2, point1, point2, xlo, xhi, ylo, yhi, adjust=2.0):
            # cut
            sheet1.cut(point1[0], point1[1]-adjust, point2[0], point2[1]-adjust, 1) # Keep upper part
            sheet1.cut(xlo, ylo, xlo, yhi, 0) # Keep righter part
            sheet1.cut(xhi, ylo, xhi, yhi, 1) # keep lefter part
            sheet1.cut(xlo, yhi, xhi, yhi, 0) # keep lower part

            sheet2.cut(point1[0], point1[1]+adjust, point2[0], point2[1]+adjust, 0) # Keep lower part
            sheet2.cut(xlo, ylo, xlo, yhi, 0) # Keep righter part
            sheet2.cut(xhi, ylo, xhi, yhi, 1) # keep lefter part
            sheet2.cut(xlo, -yhi, xhi, -yhi, 1) # keep higher part

            #sheet2.translate(0, -0.3*yy)
            sheet3 = copy.deepcopy(sheet1)
            sheet3.merge(sheet2)

            return sheet3

        def intersection(point1, point2, point3, point4):
            # Calculate the intersction of (point1, point2) and (point3, point4)
            x1 = point1[0]; y1 = point1[1]
            x2 = point2[0]; y2 = point2[1]
            x3 = point3[0]; y3 = point3[1]
            x4 = point4[0]; y4 = point4[1]
            coeff1 = np.array([[y2-y1, x1-x2], [y4-y3, x3-x4]])
            coeff2 = np.array([y2*x1-y1*x2, y4*x3-y3*x4])
            intsec = np.linalg.solve(coeff1, coeff2)
            return intsec

        theta1 = 15; theta2 = -15

        if method == 1:
            # sheet1 (up)
            sheet1 = self
            sheet1.rectangle_sheet(100, 100)
            sheet1.translate(-50*xx, -50*yy)
            sheet1.rotate(theta1)
            # sheet2 (down)
            sheet2 = Material2D()
            sheet2.copy(sheet1)
            sheet2.rotate(theta2)
            # cut and merge
            sheet3 = sheet_cut_merge(sheet1, sheet2, [-1, 0], [30*xx, 0], 0, 10*xx, 0, 50*yy, 0)
            self.atomList = sheet3.atomList

        if method == 2: 
            Theta = (theta1 - theta2)/2
            k1 = math.tan(deg2rad(30-Theta))
            k2 = -math.tan(deg2rad(30+Theta))
            dx = 6*xx
            for i in range(0, 3):
                xlo = dx * i; xhi = dx * (i+1)
                point1 = [xlo, 0]; point3 = [xhi, 0]
                point11 = [xlo + 10, k1 * 10]; point33 = [xhi + 10, k2 * 10]
                point2 = intersection(point1, point11, point3, point33)
                print(point1, point2, point3)
                # sheet1 (up)
                sheet1 = Material2D()
                sheet1.rectangle_sheet(100, 100)
                sheet1.translate(-50*xx, -50*yy)
                sheet1.rotate(theta1)
                # sheet2 (down)
                sheet2 = Material2D()
                sheet2.copy(sheet1)
                sheet2.rotate(theta2)

                sheet11 = copy.deepcopy(sheet1)
                sheet22 = copy.deepcopy(sheet2)

                sheet3 = sheet_cut_merge(sheet1, sheet2, point1, point2, xlo, point2[0], 0, 10*yy)
                self.merge(sheet3)
                sheet4 = sheet_cut_merge(sheet11, sheet22, point2, point3, point2[0], xhi, 0, 10*yy)
                self.merge(sheet4)

# -------------------------- User-defined method End --------------------------------------- #  
                


    def plot(self):
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        xx = []; yy = []; zz = []
        for atom in self.atomList:
            xx.append(atom.x); yy.append(atom.y); zz.append(atom.z)
        ax.scatter(xx, yy, zz, '3r')
        plt.xlabel('x')
        plt.ylabel('y')
        #plt.zlabel('z')  
        plt.show()


    def output(self, fname):
        if self.outputStyle == "LAMMPS":
            with open(fname, 'w') as f:
                f.write("Configuration file for 2 MoS2 Particle\n\n")
                f.write("%s atoms\n"%len(self.atomList))
                f.write("%s atom types\n"%self.ntype)
                f.write("%-12.8f %-12.8f xlo xhi\n"%(self.box[0], self.box[1]))
                f.write("%-12.8f %-12.8f ylo yhi\n"%(self.box[2], self.box[3]))
                f.write("%-12.8f %-12.8f zlo zhi\n"%(self.box[4], self.box[5]))
                f.write("\n")
                f.write("Atoms\n")
                f.write("position for MoS2 Atoms\n")
                iatom = 1
                for atom in self.atomList:
                    f.write("%d %d %d %16.8f %16.8f %16.8f\n"%(iatom, atom.type, 0, atom.x, atom.y, atom.z))
                    iatom += 1
                print("Successfully generate LAMMPS input configuration file %s"%fname)

        if self.outputStyle == "VESTA":
            with open(fname, 'w') as f:
                f.write("%s\n"%len(self.atomList))
                f.write("\n")
                f.write("%-12.8f %-12.8f xlo xhi\n"%(self.box[0], self.box[1]))
                f.write("%-12.8f %-12.8f ylo yhi\n"%(self.box[2], self.box[3]))
                f.write("%-12.8f %-12.8f zlo zhi\n"%(self.box[4], self.box[5]))
                f.write("\n")
                for atom in self.atomList:
                    if atom.type == 1:
                        f.write("S")
                    if atom.type == 2:
                        f.write("Mo")
                    f.write("%16.8f %16.8f %16.8f\n"%(atom.x, atom.y, atom.z))

                print("Successfully generate VESTA input configuration file %s"%fname)

        return

if __name__ == '__main__':
    inFileName = "myinput.in"
    # Read input file
    myinput = InputParams(inFileName)

    # Generate configurations
    # Li2CO3
    #------------------------------------------
    # Li2CO3 = Material2D()
    # Li2CO3.readFromXYZ('mp-3054_Li2CO3.xyz')
    # Li2CO3.triclincLattice(5, 5, 5)
    # Li2CO3.output_triclinc("Li2CO3.xyz")
    #------------------------------------------  
    # Diamond nanopillar
    #------------------------------------------
    # diamond = Material2D()
    # diamond.readFromXYZ('AMS_DATA.xyz')
    # r = 10
    # h = 40
    # diamond.triclincLattice(4*r, 4*r, 3*h)
    # lx = 3.5567173
    # def myCylinder(atom, norm):
    #      return diamond.cylinder_cut(2*r*lx, 2*r*lx, r*lx, 0, 2*h*lx, atom, norm) or \
    #       diamond.semi_sphere_cut(2*r*lx, 2*r*lx, 2*h*lx, r*lx, atom, norm)
    # diamond.cut3D(myCylinder, 0)
    # diamond.output_triclinc("Diamond_ref.xyz", "Diamond", isTriclinc=False)
    # strain = 0.14
    # diamond.pure_bending_mapping(r*lx/strain, 500)
    # diamond.output_triclinc("Diamond.xyz", "Diamond", isTriclinc=False)
    #------------------------------------------
    # Diamond cubic
    #------------------------------------------
    # lx = 3.5567173
    # diamond.triclincLattice(40, 40, 40)
    # diamond.translate(-20*lx, -20*lx, -20*lx)
    # diamond.rotate(45)   
    # def mySurface(atom, norm):
    #     a = 5*np.sqrt(2)
    #     return diamond.surface_cut(0, a*lx, 0, a*lx, 0, a*lx, atom, norm)
    # diamond.cut3D(mySurface, 0)
    # print "box size: %s"%(5*np.sqrt(2)*lx)
    # diamond.output_triclinc("Diamond.xyz", "Diamond", isTriclinc=False)
    #------------------------------------------  
    # Li nanopillar
    #------------------------------------------
    Li = Material2D()
    Li.readFromXYZ('./Li/Li_bcc_unit_cell.xyz')
    r = 3
    h = 15
    Li.triclincLattice(4*r, 4*r, 3*h)
    lx = 3.51
    def myCylinder(atom, norm):
         return Li.cylinder_cut(2*r*lx, 2*r*lx, r*lx, 0, 2*h*lx, atom, norm)
    Li.cut3D(myCylinder, 0)
    Li.output_triclinc("./Li/Li_nanopillar.xyz", "Li nanopillar", isTriclinc=False, charge=True)
    # strain = 0.14
    # Li.pure_bending_mapping(r*lx/strain, 500)

