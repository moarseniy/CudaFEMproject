import re
import os

#3D
class FileReader:
    def __init__(self, filename, mesh_dir, raw_mesh_dir,
                 new_dir='new_mesh', dim=3):
        # Имя файла который был создан фидесисом,
        # где хранится сетка.
        self.filename = filename
        self.mesh_dir = mesh_dir
        self.raw_mesh_dir = raw_mesh_dir
        self.new_dir = new_dir
        self.dim = dim

    def parse_nodes(self):
        # Требуется распарсить и записать Nodes
        filename = self.raw_mesh_dir + '/' + self.filename
        nodes = {}

        print('Parsing file for nodes start...')

        with open(filename, 'r') as mesh:
            for line in mesh:
                # Ищем строчку с которой начинаются Nodes,
                # сохраняем все следующие строчки, пока
                # не встретим '$' - символ начала следующего блока.
                if line.strip() == '*NODE':
                    line = mesh.readline()[1:]  # Т.к. там вначале стоит символ $.
                    line = mesh.readline()  # Чтобы сразу считывать числа.
                    while line.strip() != '$':
                        line = line.strip().split(' ')
                        line = [x for x in line if x != '']
                        nodes[line[0]] = line[1:(self.dim+1)]
                        line = mesh.readline()
                        
        #print(nodes)
        print('Parsing file for nodes end.')

        self.raw_nodes = nodes
        print('OK')

    def prepare_nodes(self):
        # Требуется преобразовать массив raw_nodes

        try:
            nodes = self.raw_nodes
        except Exception as e:
            print('Use method parse_nodes first!!')
            return

        print('Prepare and writing nodes in new file start...')

        with open(self.dir_name + '/nodes.txt', 'w') as write_file:
            write_file.write(str(len(nodes)) + '\n')
            #print('Attention dimension of mesh == {0}'.format(self.dim))
            for node in range(len(nodes.keys())):
                #k = [str(node)] + list(nodes[list(nodes.keys())[node]])
                k = list(nodes[list(nodes.keys())[node]])
                line = ' '.join(k[:self.dim+1])
                write_file.write(line + '\n')

        print('Prepare and writing nodes in new file end.')
        print('OK')

    def parse_elements(self):
        # Требуется распарсить и записать Elements

        filename = self.raw_mesh_dir + '/' + self.filename
        elements = []

        print('Parsing file for elements start...')

        with open(filename, 'r') as mesh:
            for line in mesh:

                # Ищем строчку с которой начинаются Elements,
                # сохраняем все следующие строчки, пока
                # не встретим '*END' - символ начала следующего блока.
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
            #print('Attention dimension of mesh == {0}'.format(self.dim))
            for element in elements:
                line = ' '.join(re.split('\s+', element))
                line = line.split(' ')
                for i in range(2, len(line)):
                    line[i] = str(list(nodes.keys()).index(line[i]))
                #line = ' '.join(list(str(int(line[0])-1))+line[2:-1])
                line = ' '.join(line[2:6])
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
                        #set_nodes[key] = line[0:]
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

        #print(set_nodes)                
        print('Parsing file for constraints and sets end.')

        self.raw_constraints = constraints
        self.raw_set_nodes = set_nodes
        print('OK')

    def prepare_constraints(self):
        try:
            constraints = self.raw_constraints
            set_nodes = self.raw_set_nodes
        except Exception as e:
            print('Problem!!')
            return

        print('Prepare and writing constraints and sets in new files start...')
        #count the number of constraints
        count = 0
        for constraint in constraints:
            for i in set_nodes.keys():
                if (constraint[0] == i):
                    for j in range(0, len(set_nodes[i])):
                        count+=1


        with open(self.dir_name + '/constraints.txt', 'w') as write_file:
            write_file.write(str(count) + '\n')
            for constraint in constraints:
                for i in set_nodes.keys():
                    if (constraint[0] == i):
                        for j in range(0, len(set_nodes[i])):
                            line = str(set_nodes[i][j])
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
                
        #print(constraints)
        #print(set_nodes)
                
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
                        
        #print(forces)
        #print(loads)
                     
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
        #count the number of constraints
        count = 0
        for load in loads:
            for i in set_nodes.keys():
                if load[0] == i:
                    for j in range(0, len(set_nodes[i])):
                        count += 1

        with open(self.dir_name + '/loads.txt', 'w') as write_file:
            write_file.write(str(count) + '\n')
            for load in loads:
                for i in set_nodes.keys():
                    if load[0] == i:
                        for j in range(0, len(set_nodes[i])):
                            line = str(set_nodes[i][j]) + ' '
                            for k in forces.keys():
                                if load[2] == k:
                                    #line += ' '.join(forces[k])
                                    if load[1] == '1':
                                        line += forces[k][1]
                                        line += '.0 0.0 0.0'
                                    if load[1] == '2':
                                        line += '0.0 '
                                        line += forces[k][1]
                                        line += '.0 0.0'
                                    if load[1] == '3':
                                        line += '0.0 0.0 '
                                        line += forces[k][1]
                                        line += '.0'
                            write_file.write(line + '\n')
                                    
        #print(set_nodes)
        #print(loads)
        #print(forces)
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


    def make_directory(self):
        # Делаю новую директорию, куда
        # буду сбрасывать данные связанные с этой сеткой.
        # Требуется начинать работу с этой команды.
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
        # Переносит исходный файл в новую директорию.
        # Требуется заканчивать работу этой командой.

        print('Move base file in new directory.')
        try:
            os.replace(self.raw_mesh_dir + '/' + self.filename,
                       self.dir_name + '/' + self.filename)
        except Exception as e:
            print('Some troubles with moving ' + self.filename
                  + ' to ' + self.new_dir)
