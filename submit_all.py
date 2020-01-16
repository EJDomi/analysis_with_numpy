import argparse as arg
import os

if __name__ == "__main__":
     
    parser = arg.ArgumentParser()
    parser.add_argument('-d', '--submit_dir', type=str)

    args = parser.parse_args()

    submit_dir = args.submit_dir

    submit_list = [os.path.join(submit_dir, f) for f in os.listdir(submit_dir) if (os.path.isfile(os.path.join(submit_dir, f)) and ('.submit' in f))]

    for f in submit_list:
        os.system('condor_submit ' + f)

