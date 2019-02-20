import time
import subprocess
def batchShellTask(all_task,limit_task):
    '''
    batch shell task
    :param all_task: all shell task list
    :param limit_task: limit task number simultaneously running
    :type all_task: list
    :type limit_task: int
    :return: null
    '''
    task_pool=[]
    task_remain=len(all_task)
    for task in all_task:
        task_remain+=-1
        break_out = True
        p = subprocess.Popen(task, shell=True, stdin=subprocess.PIPE,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
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
def test_demo():
    shells=['ls','pwd','ifconfig']
    batchShellTask(shells,2)

if __name__=='__main__':
    test_demo()
