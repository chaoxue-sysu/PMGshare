import time
import subprocess
def batchShellTask(all_task,limit_task,Log_file):
    '''
    batch shell task
    :param all_task: all shell task list
    :param limit_task: limit task number simultaneously running
    :type all_task: list
    :type limit_task: int
    :return: null
    '''
    log = open(Log_file, "w")
    task_pool=[]
    task_remain=len(all_task)
    for task in all_task:
        task_remain+=-1
        break_out = True
        p = subprocess.Popen(task, shell=True, stdin=log, stdout=log, stderr=log, close_fds=True)
        task_pool.append(p)
        print(time.strftime("[%Y-%m-%d %H:%M:%S]", time.localtime()) + ' '+str(p.pid)+': '+task+' start ...')
        if len(task_pool)==limit_task or task_remain==0:
            while break_out:
                for intask_Popen in task_pool:
                    if intask_Popen.poll()!=None:
                        print(time.strftime("[%Y-%m-%d %H:%M:%S]", time.localtime()) + ' '+str(intask_Popen.pid)+': '+' finish...')
                        task_pool.remove(intask_Popen)
                        break_out = False
                        if task_remain==0:
                            break_out=True
                        if len(task_pool)==0:
                            break_out=False
                        break
    log.close()

def test_demo():
    parallel_shells=['ls','pwd','ifconfig']
    # Please specify the log file path here
    batchShellTask(parallel_shells,2,Log_file)

if __name__=='__main__':
    test_demo()
