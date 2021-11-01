from prepare_mesh_lib.file_reader import FileReader
import time

def main():
    print('-' * 25)
    print('Put full filename to parse mesh from:')
    filename = input()  # Файл, который надо распарсить. Полное имя.
    # mesh_dir = input()
    mesh_dir = 'prepared_meshes'  # Папка где будут лежать готовые сетки.
    # raw_mesh_dir = input()
    raw_mesh_dir = 'raw_meshes'  # Папка где лежат файлы с неготовыми сетками.
    # new_dir = input()
    new_dir = filename[:-2]  # Папка куда сложить результат.


    parser = FileReader(filename, mesh_dir, raw_mesh_dir, new_dir)
    parser.make_directory()

    parser.parse_nodes()
    parser.prepare_nodes()

    parser.parse_elements()
    parser.prepare_elements()

    parser.parse_constraints_and_sets()
    parser.prepare_constraints()

    parser.parse_loads()
    parser.prepare_loads()

    #t0 = time.time()
    #parser.make_colors()
    #t1 = time.time() - t0
    #print(t1)

    # parser.make_good()

    print('-' * 25)
    print('Done')


if __name__ == '__main__':
    main()
