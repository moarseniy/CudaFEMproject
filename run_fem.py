from fem_parser import run_parser
import os, json, sys

if __name__ == '__main__':
    try:
        with open(str(sys.argv[1]), 'r') as file_read:
            config_data = json.load(file_read)
        run_parser(str(config_data['task']),int(config_data['DIM']))
    except Exception as e:
        print('Something wrong with config file')
        sys.exit(0)

    
