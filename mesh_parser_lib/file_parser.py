import re, os

class FileParser:
    def __init__(self, filename, mesh_dir, raw_mesh_dir,
                 new_dir, dim):
        self.filename = filename
        self.mesh_dir = mesh_dir
        self.raw_mesh_dir = raw_mesh_dir
        self.new_dir = new_dir
        self.dim = dim

    def make_directory(self):
        print('Make new directory named: ' + self.new_dir
              + ' in ' + self.mesh_dir)

        self.dir_name = self.mesh_dir + '/' + self.new_dir
        try:
            os.mkdir(self.dir_name)
        except Exception as e:
            print('This directory already exist!')
            return
        print('OK')

    def make_good(self):
        print('Move base file in new directory.')
        try:
            os.replace(self.raw_mesh_dir + '/' + self.filename,
                       self.dir_name + '/' + self.filename)
        except Exception as e:
            print('Some troubles with moving ' + self.filename
                  + ' to ' + self.new_dir)
            
    def parse_nodes(self):
        filename = self.raw_mesh_dir + '/' + self.filename
        nodes = {}

        print('Parsing file for nodes start...')

        with open(filename, 'r') as mesh:
            for line in mesh:
                if line.strip() == '*NODE':
                    line = mesh.readline()                
                    line = mesh.readline()
                    while line.strip() != '$':
                        line = line.strip().split(' ')
                        line = [x for x in line if x != '']
                        nodes[line[0]] = line[1:(self.dim+1)]
                        line = mesh.readline()
        print('Parsing file for nodes end.')
        self.raw_nodes = nodes
        print('OK')

    def prepare_nodes(self):
        try:
            nodes = self.raw_nodes
        except Exception as e:
            print('Use method parse_nodes first!!')
            return

        print('Prepare and writing nodes in new file start...')
        with open(self.dir_name + '/nodes.txt', 'w') as write_file:
            write_file.write(str(len(nodes)) + '\n')
            for node in range(len(nodes.keys())):
                k = list(nodes[list(nodes.keys())[node]])
                line = ' '.join(k[:self.dim+1])
                write_file.write(line + '\n')

        print('Prepare and writing nodes in new file end.')
        print('OK')

    def parse_elements(self):
        filename = self.raw_mesh_dir + '/' + self.filename
        elements = []

        print('Parsing file for elements start...')

        with open(filename, 'r') as mesh:
            for line in mesh:
                if line.strip() == '*ELEMENT_SOLID':
                    line = mesh.readline()
                    while line.strip()[0] != '*':
                        elements.append(line.strip())
                        line = mesh.readline()

        print('Parsing file for elements end.')
        self.raw_elements = elements
        print('OK')

    def prepare_elements(self):
        try:
            elements = self.raw_elements
            nodes = self.raw_nodes
        except Exception as e:
            print('Use mehtod parse_elements first!!')
            return

        print('Prepare and writing elements in new file start...')

        with open(self.dir_name + '/elements.txt', 'w') as write_file:
            write_file.write(str(len(elements)) + '\n')
            for element in elements:
                line = ' '.join(re.split('\s+', element))
                line = line.split(' ')[2:]
                for i in range(0, len(line)):
                    line[i] = str(list(nodes.keys()).index(line[i]))
                line = ' '.join(line[0:4])
                write_file.write(line + '\n')

        print('Prepare and writing elements in new file end.')
        print('OK')

    def parse_constraints_and_sets(self):

        filename = self.raw_mesh_dir + '/' + self.filename
        constraints = []
        set_nodes = {}
        nodes = self.raw_nodes

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
                            line[i] = str(list(nodes.keys()).index(line[i]))
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

        print('Prepare and writing constraints and sets in new files start...')
        with open(self.dir_name + '/constraints.txt', 'w') as write_file:
            for constraint in constraints:
                for set_num in set_nodes.keys():
                    if (constraint[0] == set_num):
                        write_file.write(str(len(set_nodes[set_num])) + '\n')
                        for j in range(0, len(set_nodes[set_num])):
                            line = str(set_nodes[set_num][j])
                            
                            if self.dim == 3:
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
                
        print('Prepare and writing constraints and sets in new files end.')
        print('OK')

    def parse_loads(self):
        filename = self.raw_mesh_dir + '/' + self.filename
        forces = {}
        loads = []

        print('Parsing file for loads start...')

        with open(filename, 'r') as mesh:
            for line in mesh:
                if line.strip() == '*LOAD_NODE_SET':
                    line = mesh.readline().strip().split(' ')
                    line = [x for x in line if x != '']
                    loads.append(line[:])

        with open(filename, 'r') as mesh:
            for line in mesh:
                if line.strip() == '*DEFINE_CURVE':
                    line = mesh.readline().strip()
                    key = line
                    line = mesh.readline()
                    while line.strip()[0] != '*':
                        line = line.strip().split(' ')
                        line = [x for x in line if x != '']
                        forces[key] = line[0:2]
                        forces[key].append(str(0.0))
                        line = mesh.readline()    
                     
        print('Parsing file for loads end.')
        self.raw_forces = forces
        self.raw_loads = loads
        print('OK')

    def prepare_loads(self):
        try:
            loads = self.raw_loads
            forces = self.raw_forces
            set_nodes = self.raw_set_nodes
        except Exception as e:
            print('Problem!!')
            return

        print('Prepare and writing loads in new files start...')
        with open(self.dir_name + '/loads.txt', 'w') as write_file:
            for load in loads:
                for set_num in set_nodes.keys():
                    if load[0] == set_num:
                        write_file.write(str(len(set_nodes[set_num])) + '\n')
                        for j in range(0, len(set_nodes[set_num])):
                            line = str(set_nodes[set_num][j]) + ' '
                            for force_num in forces.keys():
                                if load[2] == force_num:
                                    if load[1] == '1':
                                        line += str(forces[force_num][1]) + '.0 0.0 0.0'
                                    if load[1] == '2':
                                        line += '0.0 ' + str(forces[force_num][1]) + '.0 0.0'
                                    if load[1] == '3':
                                        line += '0.0 0.0 ' + str(forces[force_num][1]) + '.0'
                            write_file.write(line + '\n')
                                    
        print('Prepare and writing loads in new files end.')
        print('OK')

    def make_colors(self):
        raw_noads = self.raw_nodes
        raw_elements = self.raw_elements
        elements = []
        elementsColor = []
       
        for i in range(len(raw_elements)):
            temp = raw_elements[i].split(' ')
            temp = [x for x in temp if x != '']
            del temp[1]
            elements.append(temp[1:5])
        print(elements)

        k_prev, k, i = 0, 0, 0 
        size = len(elements)
        while len(elementsColor) != size:
            while i < len(elements):
                add = True
                for j in range(k_prev, k + k_prev):
                    for t in range(4):
                        for p in range(4):
                            if elements[i][t] == elementsColor[j][p] or elements[i][t] == elementsColor[j][p] or \
                            elements[i][t] == elementsColor[j][p] or elements[i][t] == elementsColor[j][p]:
                                add = False
                if add == True:
                    elementsColor.append(elements[i])
                    del elements[i]
                    k += 1
                else:
                    i += 1
            print(k)
            k_prev += k
            k = 0
            i = 0
