# -*- coding: utf-8 -*-
"""
Created on Thu Jul 21 17:27:10 2016

@author: mads

This script download all images from roboweedsupport.com and extracts the thumbnails
Data has to be exportet beforehand
"""

from concurrent.futures import ThreadPoolExecutor
import requests
import os
import sys
import multiprocessing

listOfUploadDirs=range(2,8)


Ncpus=multiprocessing.cpu_count()

def main(listOfUploadDirs):
    executor = ThreadPoolExecutor(max_workers=Ncpus)
    futures = []
    
    currentPath=os.path.split(os.path.realpath(__file__))[0]
        
    for i in listOfUploadDirs:
        if not os.path.exists('export/'+str(i)):
            #Try to download the export data, if success call the main download script
            #status = call(["wget", 'http://roboweedsupport.com/ExportedData/data_'+str(i)+'.zip'])
            #filedir='/media/mads/79131cf3-d458-45c7-bb29-830bc120a265/Phd/Software/System Software/Python/GetThumbnailsFromServer/'
            
            filedir=currentPath
            
            r = requests.head('http://roboweedsupport.com/ExportedData/data_'+str(i)+'.zip')
            if r.status_code==200: #200 = request accepted
                cmd='python "'+os.path.join(filedir,'GetThumbnailsFromServer.py')+'" --uploadid '+str(i)
                
            	#Start a new python interpreter for this upload
                a = executor.submit(os.system, cmd)
                futures.append(a)
            else:
                print('data does not exist')
            
    executor.shutdown(wait=True)


if __name__ == '__main__':
    # Arguments - integers, except for binary/continous. Default uses binary.
    # Run with option --continuous for continuous output.
    import argparse
    parser = argparse.ArgumentParser(description='Command line options')    
    #parser.set_defaults(binary=True)    
    #parser.add_argument('--uploadid', type=str, dest='listOfUploadDirs', default='1,250')
    parser.add_argument('-i','--listOfUploadDirs', nargs='+', help='<Required> Set flag', required=False, default=[2,3,4,5])
    args = parser.parse_args(sys.argv[1:])
    main(**{k:v for (k,v) in vars(args).items() if v is not None})
