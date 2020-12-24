import threading
import os
from concurrent.futures import ThreadPoolExecutor
import time, datetime

def test(value1, value2=None):
    print("%s threading is printed %s, %s"%(threading.current_thread().name, value1, value2))
#     time.sleep(2)

def run_obabel(src, dst):
    cmd = 'obabel ' + str(src) + ' -O ' + dst
    os.system(cmd)

def run_obabel_multithread(input_dir, output_dir, max_workers = 100):
    start_time = datetime.datetime.now()
    print("start time: " + str(start_time))

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    input_smis = []
    output_paths = []
    if os.path.isdir(input_dir):
        for file in os.listdir(input_dir):
            file_path = os.path.join(input_dir, file)
            ouput_path = os.path.join(output_dir, file + '.png')
            output_paths.append(ouput_path)
            with open(file_path, 'r', encoding='utf8') as f:
                smi = f.readlines()[0].replace('\n', '')
                input_smis.append(smi)
    else:
        print('input error')

    threadPool = ThreadPoolExecutor(max_workers=max_workers, thread_name_prefix="execute_obabel_")
    threadPool.map(run_obabel, input_smis, output_paths)
    threadPool.shutdown(wait=True)

    end_time = datetime.datetime.now()
    print("end time: " + str(datetime.datetime.now()))
    print('cost time:' + str(end_time - starttime))

if __name__ == "__main__":
    starttime = datetime.datetime.now()
    threadPool = ThreadPoolExecutor(max_workers=4, thread_name_prefix="test_")
    value1 = [i for i in range(1000)]
    value2 = [i+1 for i in range(1000)]
    # for i in range(0,100):
    #     threadPool.map(test, [i],[i+1])
    threadPool.map(test, value1, value2)
    threadPool.shutdown(wait=True)
    endtime = datetime.datetime.now()
    print('cost time:' + str(endtime - starttime))