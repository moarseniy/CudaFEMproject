import re, os, time
import math
import numpy as np
import collections

# Write a line to the beginning of file. Used in prepare_constraints() method
def line_prepender(filename, line):
    with open(filename, 'r+') as f:
        content = f.read()
        f.seek(0, 0)
        f.write(line.rstrip('\r\n') + '\n' + content)

class FileParser:
    def __init__(self, filename, prepared_mesh_dir, raw_mesh_dir,
                 task_name, dim):
        self.filename = filename
        self.prepared_mesh_dir = prepared_mesh_dir
        self.raw_mesh_dir = raw_mesh_dir
        self.task_name = task_name
        self.dim = dim

    def make_directory(self):
        print('Make new directory named: ' + self.task_name
              + ' in ' + self.prepared_mesh_dir)

        self.prepared_mesh_dir = self.prepared_mesh_dir
        try:
            os.mkdir(self.prepared_mesh_dir)
        except Exception as e:
            print('This directory already exists!')
            return
        print('OK')

    def parse_nodes(self):
        filename = self.raw_mesh_dir
        nodes = []
        nodes_num = {}

        print('Parsing file for nodes start...')

        with open(filename, 'r') as mesh:
            for line in mesh:
                if line.strip() == '*NODE':
                    line = mesh.readline()                
                    line = mesh.readline()
                    while line.strip() != '$':
                        line = line.strip().split(' ')
                        line = [x for x in line if x != '']
                        nodes.append(line[:(self.dim+1)])
                        line = mesh.readline()
        nodes.sort(key=lambda x: int(x[0]))
        for i in range(len(nodes)):
            nodes_num[nodes[i][0]] = i
            #nodes[i] = nodes[i][1:]
            
        print('Parsing file for nodes end.')
        self.nodes_num = nodes_num
        self.raw_nodes = nodes
        print('OK')

    def prepare_nodes(self):
        try:
            nodes = self.raw_nodes
            dim = self.dim
        except Exception as e:
            print('Use method parse_nodes first!!')
            return

        print('Prepare and writing nodes in new file start...')
        with open(self.prepared_mesh_dir + '/nodes.txt', 'w') as write_file:
            write_file.write(str(len(nodes)) + '\n')
            for el in nodes:
                line = ' '.join(el[1:dim+1])
                write_file.write(line + '\n')

        print('Prepare and writing nodes in new file end.')
        print('OK')

    def parse_elements(self):
        filename = self.raw_mesh_dir
        elements = []

        print('Parsing file for elements start...')

        with open(filename, 'r') as mesh:
            for line in mesh:
                if line.strip() == '*ELEMENT_SOLID' or line.strip() == '*ELEMENT_SHELL':
                    line = mesh.readline()
                    while line.strip()[0] != '*':
                        elements.append(' '.join(re.split('\s+', line.strip())).split(' '))
                        line = mesh.readline()

        print('Parsing file for elements end.')
        self.raw_elements = elements
        print('OK')

    def prepare_elements(self):
        try:
            elements = self.raw_elements
            nodes = self.raw_nodes
            nodes_num = self.nodes_num
            dim = self.dim
            load_segments = self.raw_load_segments
        except Exception as e:
            print('Use mehtod parse_elements first!!')
            return

        print('Prepare and writing elements in new file start...')
        start_time = time.time()

        with open(self.prepared_mesh_dir + '/elements.txt', 'w') as write_file:
            write_file.write(str(len(elements)) + '\n')
            for element in elements:
                line = element[2:]
                line = line[:dim+1]
                if dim == 2:
                    for i in range(len(load_segments)):
                        if load_segments[i][3] in line and load_segments[i][4] in line:
                            idx = 3 - line.index(load_segments[i][3]) - line.index(load_segments[i][4])
                            # int((len(line)-1)*len(line)/2) = 3
                            load_segments[i].append(line[idx])
                            load_segments[i].append(elements.index(element))
                if dim ==3:
                    for i in range(len(load_segments)):
                        if load_segments[i][3] in line and load_segments[i][4] in line and load_segments[i][5] in line:
                            idx = 6 - line.index(load_segments[i][3]) - line.index(load_segments[i][4]) - line.index(load_segments[i][5])
                            # int((len(line)-1)*len(line)/2) = 6
                            load_segments[i].append(line[idx])
                            load_segments[i].append(elements.index(element))
                for i in range(len(line)):
                    line[i] = str(nodes_num[line[i]])
                line = ' '.join(line[:])
                write_file.write(line + '\n')
        
        print('Prepare and writing elements in new file end. ' + str(time.time() - start_time))
        self.raw_load_segments = load_segments
        print('OK')

    def parse_constraints_and_sets(self):

        filename = self.raw_mesh_dir
        constraints = []
        set_nodes = {}
        nodes_num = self.nodes_num

        print('Parsing file for constraints and sets start...')

        with open(filename, 'r') as mesh:
            for line in mesh:            
                if line.strip() == '*SET_NODE_LIST':
                    line = mesh.readline().strip()
                    key = line
                    line = mesh.readline()
                    temp = []
                    while line.strip()[0] != '*':
                        line = line.strip().split(' ')
                        line = [x for x in line if x != '']
                        for i in range(0, len(line)):
                            line[i] = str(nodes_num[line[i]])
                            temp.append(line[i])
                        set_nodes[key] = temp
                        line = mesh.readline()

        with open(filename, 'r') as mesh:
            for line in mesh:
                if line.strip() == '*BOUNDARY_SPC_SET':
                    line = mesh.readline().strip().split(' ')
                    line = [x for x in line if x != '']
                    constraints.append(line[:])
             
        print('Parsing file for constraints and sets end.')
        self.raw_constraints = constraints
        self.raw_set_nodes = set_nodes
        print('OK')

    def prepare_constraints(self):
        try:
            constraints = self.raw_constraints
            set_nodes = self.raw_set_nodes
            dim = self.dim
        except Exception as e:
            print('Problem!!')
            return

        num_constraints = 0
        print('Prepare and writing constraints and sets in new files start...')
        with open(self.prepared_mesh_dir + '/constraints.txt', 'w') as write_file:
            for constraint in constraints:
                for set_num in set_nodes.keys():
                    if (constraint[0] == set_num):
                        num_constraints += len(set_nodes[set_num])
                        for j in range(0, len(set_nodes[set_num])):
                            line = str(set_nodes[set_num][j])
                            if dim == 2:
                                if int(constraint[2])==1 and int(constraint[3])==0:
                                    line += ' ' + str(1) + '\n'
                                elif int(constraint[2])==0 and int(constraint[3])==1:
                                    line += ' ' + str(2) + '\n'
                                elif int(constraint[2])==1 and int(constraint[3])==1:
                                    line += ' ' + str(3) + '\n'
                            if dim == 3:
                                if int(constraint[2])==1 and int(constraint[3])==0 and int(constraint[4])==0:
                                    line += ' ' + str(1) + '\n'
                                elif int(constraint[2])==0 and int(constraint[3])==1 and int(constraint[4])==0:
                                   line += ' '+ str(2) + '\n'
                                elif int(constraint[2])==0 and int(constraint[3])==0 and int(constraint[4])==1:
                                    line += ' ' + str(3) + '\n'
                                elif int(constraint[2])==1 and int(constraint[3])==1 and int(constraint[4])==0:
                                    line += ' ' + str(4) + '\n'
                                elif int(constraint[2])==1 and int(constraint[3])==0 and int(constraint[4])==1:
                                    line += ' ' + str(5) + '\n'
                                elif int(constraint[2])==0 and int(constraint[3])==1 and int(constraint[4])==1:
                                    line += ' ' + str(6) + '\n'
                                elif int(constraint[2])==1 and int(constraint[3])==1 and int(constraint[4])==1:
                                    line += ' ' + str(7) + '\n'
                            write_file.write(line)

        line_prepender(self.prepared_mesh_dir + '/constraints.txt', str(num_constraints))        
        print('Prepare and writing constraints and sets in new files end.')
        print('OK')

    def parse_loads(self):
        filename = self.raw_mesh_dir
        forces = {}
        load_node_sets = []
        load_segments = []

        print('Parsing file for loads start...')

        with open(filename, 'r') as mesh:
            for line in mesh:
                if line.strip() == '*LOAD_NODE_SET':
                    line = mesh.readline().strip().split(' ')
                    line = [x for x in line if x != '']
                    load_node_sets.append(line[:])
                    
        with open(filename, 'r') as mesh:
            for line in mesh:
                if line.strip() == '*LOAD_SEGMENT':
                    line = mesh.readline().strip().split(' ')
                    line = [x for x in line if x != '']
                    load_segments.append(line[:])
        
        with open(filename, 'r') as mesh:
            for line in mesh:
                if line.strip() == '*DEFINE_CURVE':
                    line = mesh.readline().strip()
                    key = line
                    line = mesh.readline()
                    while line.strip()[0] != '*':
                        line = line.strip().split(' ')
                        line = [x for x in line if x != '']
                        forces[key] = line[1]
                        line = mesh.readline()
                        
        print('Parsing file for loads end.')
        self.raw_forces = forces
        self.raw_load_node_sets = load_node_sets
        self.raw_load_segments = load_segments
        print('OK')

    def prepare_loads(self):
        try:
            load_node_sets = self.raw_load_node_sets
            load_segments = self.raw_load_segments
            forces = self.raw_forces
            set_nodes = self.raw_set_nodes
            dim = self.dim
            nodes = self.raw_nodes
            nodes_num = self.nodes_num
        except Exception as e:
            print('Problem in loads!!')
            return
        print('Prepare and writing loads in new files start...')
        
        with open(self.prepared_mesh_dir + '/loads.txt', 'w') as write_file:
            if len(load_node_sets) == 0:
                write_file.write(str(0))
            for load in load_node_sets:
                for set_num in set_nodes.keys():
                    if load[0] == set_num:
                        write_file.write(str(len(set_nodes[set_num])) + '\n')
                        for j in range(len(set_nodes[set_num])):
                            line = str(set_nodes[set_num][j]) + ' '
                            for force_num in forces.keys():
                                if load[2] == force_num:
                                    if dim == 2:
                                        if load[1] == '1':
                                            line += str(float(forces[force_num])) + ' 0.0'
                                        if load[1] == '2':
                                            line += '0.0 ' + str(float(forces[force_num]))
                                    elif dim == 3:
                                        if load[1] == '1':
                                            line += str(float(forces[force_num])) + ' 0.0 0.0'
                                        if load[1] == '2':
                                            line += '0.0 ' + str(float(forces[force_num])) + ' 0.0'
                                        if load[1] == '3':
                                            line += '0.0 0.0 ' + str(float(forces[force_num]))
                            write_file.write(line + '\n')
                            
        with open(self.prepared_mesh_dir + '/stress.txt', 'w') as write_file:
            write_file.write(str(len(load_segments)) + '\n')
            for load in load_segments:
                if dim == 2:
                    x1 = float(nodes[nodes_num[load[3]]][0])
                    y1 = float(nodes[nodes_num[load[3]]][1])
                    
                    x2 = float(nodes[nodes_num[load[4]]][0])
                    y2 = float(nodes[nodes_num[load[4]]][1])
                    
                    xb = float(nodes[nodes_num[load[8]]][0])
                    yb = float(nodes[nodes_num[load[8]]][1])

                    xn = y2 - y1  #/ (math.sqrt(y1*y1+x1*x1))
                    yn = -(x2 - x1) #/ (math.sqrt(y1*y1+x1*x1))

                    dot_prod = xn * (xb - x1) + yn * (yb - y1)
                    if dot_prod > 0:
                        xn = -xn
                        yn = -yn
                    
                    for force_num in forces.keys():
                        if load[0] == force_num:
                            pressure = float(forces[force_num])
                            line = str(nodes_num[load[3]]) + ' ' + \
                                   str(nodes_num[load[4]]) + ' ' + \
                                   str(nodes_num[load[8]]) + ' ' + \
                                   str(load[9]) + ' ' + str(xn) + ' ' + str(yn) + ' ' + str(pressure) + '\n'
                            #line += str(pressure * xn1) + ' ' + str(pressure * yn1) + '\n'
                            #line += str(list(nodes.keys()).index(load[4])) + ' '
                            #line += str(pressure * xn2) + ' ' + str(pressure * yn2) + '\n'
                            write_file.write(line)
                    # TODO: DISCUSS HOw TO IMPROVE DIFFERENT CASES
                elif dim == 3:
                    
                    x1 = float(nodes[nodes_num[load[3]]][0])
                    y1 = float(nodes[nodes_num[load[3]]][1])
                    z1 = float(nodes[nodes_num[load[3]]][2])
                    
                    x2 = float(nodes[nodes_num[load[4]]][0])
                    y2 = float(nodes[nodes_num[load[4]]][1])
                    z2 = float(nodes[nodes_num[load[4]]][2])

                    x3 = float(nodes[nodes_num[load[5]]][0])
                    y3 = float(nodes[nodes_num[load[5]]][1])
                    z3 = float(nodes[nodes_num[load[5]]][2])

                    xb = float(nodes[nodes_num[load[8]]][0])
                    yb = float(nodes[nodes_num[load[8]]][1])
                    zb = float(nodes[nodes_num[load[8]]][2])

                    a = [x2 - x1, y2 - y1, z2 - z1]
                    b = [x3 - x1, y3 - y1, z3 - z1]

                    n = np.cross(a, b)

                    dot_prod = n[0] * (xb - x1) + n[1] * (yb - y1) + n[2] * (zb - z1)
                    if dot_prod > 0:
                        n[0] = -n[0]
                        n[1] = -n[1]
                        n[2] = -n[2]

                    for force_num in forces.keys():
                        if load[0] == force_num:
                            pressure = float(forces[force_num])
                            line = str(nodes_num[load[3]]) + ' ' + \
                                   str(nodes_num[load[4]]) + ' ' + \
                                   str(nodes_num[load[5]]) + ' ' + \
                                   str(load[9]) + ' ' + str(n[0]) + ' ' + str(n[1]) + \
                                   ' ' + str(n[2]) + ' ' + str(pressure) + '\n'
                            write_file.write(line)
            
        print('Prepare and writing loads in new files end.')
        print('OK')
