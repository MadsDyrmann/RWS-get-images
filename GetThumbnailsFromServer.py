# -*- coding: utf-8 -*-
"""
Created on Wed Dec 16 10:12:00 2015

@author: mads
"""
from __future__ import print_function
import sys
import os
import urllib2
import zipfile
#import scipy.ndimage
#import scipy.misc
import json
import cv2
import numpy as np
import csv
import re
import codecs
import time
import PIL
#from PIL import PngImagePlugin
from PIL import Image
import PIL.ExifTags
from QRreader.QRArea import QRreader
#from selenium import webdriver
from threading import Thread
jsonreader=json.JSONDecoder()
import requests
from concurrent.futures import ThreadPoolExecutor
from csv import reader


def main(uploadid='0',skipUnsorted=False, exportThumbs=True,exportKitti = True, exportYOLO = False, exportFastRCNN = False):

    print('Processing id no. '+uploadid)    
    
    #Get the annotation data
    getRoboWeedSupportAnnotationData(uploadid)    

    #See if there are any annotations
    anyValidAnnotation=anyAnnotations(uploadid)

    #If there are at least one annotation, proceed with downloading the images
    if anyValidAnnotation:
        
        ##get the images based on valid annotations
        #getImagesBasedOnAnnotations(uploadid,allAnnotations=True)
        
        ##Get all images
        #getRoboWeedSupportUpload(uploadid)
        
        #Use the API to get list of images
        getImagesBasedOnImageAPI(uploadid)

        #get list of images
        onlyFiles = getImagelist(uploadid)
    
        #get list of species in PVO
        speciesDictionary = getSpeciesDirectory()
    
        for i in range(len(onlyFiles)):
            print('Image',i+1,'of',len(onlyFiles),'in upload id',uploadid)
            #find the annotation, if it exist
            (dirname, filename) = os.path.split(onlyFiles[i])
            filename, file_extension = os.path.splitext(filename)
            annotationfile=os.path.join(dirname,'data_'+filename+'.txt')
            if os.path.exists(annotationfile):
    
                #Read the annotation
                with open(annotationfile, 'r') as f:
                    annotations = f.readlines()
                
                #Convert to dictionary of annotations
                annotationDict = getAnnotationDictionary(annotations, speciesDictionary)
                
                #Remove unsorted and Junk
                if skipUnsorted:
                    removeplants=['Unsorted','Junk']
                    for label in removeplants:
                        annotationDict=[anno for anno in annotationDict if label not in anno['species']]
    
                #Export thumbnails
                if exportThumbs:
                    exportThumbnails(annotationDict, dirname, onlyFiles[i])
                if exportFastRCNN:
                    exportLocalizationDataforFastRCNN(imagepath=onlyFiles[i],annotations=annotationDict)
                if exportYOLO:
                    exportLocalizationDataForYOLO(imagepath=onlyFiles[i],annotations=annotationDict)
                if exportKitti:
                    exportDetectionDataKittiFormat(imagepath=onlyFiles[i],annotations=annotationDict)

#This method checks if the are any annotations for a given upload
def anyAnnotations(uploadid):
    annotationfilelist = []
    for root, dirnames, filenames in os.walk(uploadid):
        for filename in filenames:
            if filename.lower().endswith(('.txt')):
                annotationfilelist.append(os.path.join(root, filename))    
    
    for annotationFile in annotationfilelist:
        #Read the annotation
        with open(annotationFile, 'r') as f:
            annotations = f.readlines()
            
            if len(annotations):
                return True
    return False #Only reached if there is no valid annotations



def getSpeciesDirectory():
    speciesDictionary={}
    currentDir=os.path.dirname(os.path.realpath(__file__))
    ukrudtdir=os.path.join(currentDir,'Ukrudtsarter i PVO til IGIS.CSV')
    with codecs.open(ukrudtdir,'r') as csvfile:
    #with codecs.open('Ukrudtsarter i PVO til IGIS.CSV','r') as csvfile:
        reader = csv.reader(csvfile,delimiter=';',skipinitialspace=True)
        for row in reader:
            speciesDictionary[row[0]]=row[1]    
    return speciesDictionary


def getImagelist(imagedir):
    onlyFiles = []
    for root, dirnames, filenames in os.walk(imagedir):
        for filename in filenames:
            if filename.lower().endswith(('.png','.jpg','.jpeg','.tif','.tiff','.gif')):
                onlyFiles.append(os.path.join(root, filename))
    return onlyFiles

def createIfNotExist(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)    

#Download a zip with all images form a single upload
def getRoboWeedSupportUpload(uploadid):
    #See if folder has been downloaded previously
    if not os.path.exists('./'+uploadid):
        #If folder and zip-file don't exist
        url='http://roboweedsupport.com/Download/' + uploadid + '.zip'
        if not os.path.exists('./'+uploadid + ".zip"): #if zip does not exist, download it
            os.system('wget -q '+url)
            
        #If folder doesn't exist, but zip-file does. Check if age of the zip file and see if a now should be downloaded
        else:
            f = urllib2.urlopen(url)
            fileinfo=f.headers.items()
            # NB Ved at give wget parameteren -N, sørger den automatisk for kun at downloade, hvis online-filen er nyere
            #Test for differece in file size
            onlineFileModifiedDate=[x[1] for x in fileinfo if 'last-modified' in x[0]][0]
            onlineFileModifiedDate = time.mktime(time.strptime(onlineFileModifiedDate, "%a, %d %b %Y %H:%M:%S %Z")) #time since epoch
            localFileModifiedDate=os.path.getmtime('./'+uploadid + ".zip") #time since epoch
            timedelta=onlineFileModifiedDate-localFileModifiedDate #time difference in seconds
            print('Difference in age of zip-file is ' + time.strftime('%H:%M:%S', time.gmtime(timedelta)))
            
            #Test for differences in size of the files
            onlineFileSize=[x[1] for x in fileinfo if 'content-length' in x[0]][0]
            onlineFileSize=int(onlineFileSize)
            localFilesize=os.path.getsize('./'+uploadid + ".zip")
            
            #if the file is more than 24 hours old, download it anyways
            if abs(1.0*timedelta)/3600>24:
                print ('difference in hours is '+str(abs(1.0*timedelta)/3600))
                os.system('wget -q -N '+url)
            elif onlineFileSize>localFilesize: #If size is different (download was canceled)
                print ('difference in size is',abs(onlineFileSize-localFilesize))
                os.system('wget -q -N '+url)
            

        #Unzip if successfully downloads, otherwise force download using slow hack method
        if os.path.exists('./'+uploadid + ".zip"):
            url=uploadid + ".zip"
            unzipFile(url,path=uploadid)
        else:
            getRoboWeedSupportUploadHack(uploadid)

    


def getRoboWeedSupportAnnotationData(uploadid):
     
    #Get the annotation data
    if not os.path.exists('./' + "data_"+uploadid +".zip"): #if zip does not exist, download it
        #Old method,
#        webpage = r"http://roboweedsupport.com/Administration/ViewUpload?uploadid="+str(uploadid) # edit me
#        
#        browser = webdriver.Firefox()
#        #browser.implicitly_wait(5)
#        browser.get(webpage)
#        browser.find_element_by_class_name('btn-primary').click() 
#        
#        browser.close()    
        #New method
        pass
    requests.get('http://roboweedsupport.com/api/export/annotation/'+str(uploadid))

    url='http://roboweedsupport.com/ExportedData/data_'+uploadid+'.zip'
    #wget.download(url,bar=wget.bar_adaptive)
    localfile='data_' + uploadid + ".zip"
    
    os.system('wget -N '+url+' -O "'+localfile+'"')
        

    url='data_' + uploadid + ".zip"
    unzipFile(localfile,path=uploadid)




#Download all images using API
def getImagesBasedOnImageAPI(uploadid):

    imagelistUrl='http://roboweedsupport.com/api/upload/list/'+str(uploadid)
    
    #Get list of images as string
    file = urllib2.urlopen(imagelistUrl)
    data = file.read()
    file.close()
   
    #Parse string to list of urls for images
    datalist=[line[0] for line in reader(data) if any(line) and len(line)==1 and len(line[0])>5]
    

    
    #download file
    #Download in parallel
    DownladImagesExecutor = ThreadPoolExecutor(max_workers=10)
    Downloadfutures = []    

    baseimagepath='http://roboweedsupport.com/Data' #'http://roboweedsupport.com/Data/233/20160518_135608.jpg'
    for imageurl in datalist:
        filename=os.path.split(imageurl)[-1]
        a = DownladImagesExecutor.submit(getImage, baseimagepath,uploadid,filename)
        Downloadfutures.append(a)           

    DownladImagesExecutor.shutdown(wait=True)





#Download image if there are any annotations
#This function reads the name of the annotationfiles, and if there is an annotation for a given image
#the image is downloaded. That way, only images with annotations are downloaded
# if allAnnotations is true, an image is downloaded even if there are no annotations. This is useful for getting the barcode
def getImagesBasedOnAnnotations(uploadid,allAnnotations=False):
    onlyFiles = []
    for root, dirnames, filenames in os.walk(uploadid):
        for filename in filenames:
            if filename.lower().endswith(('.txt')):
                onlyFiles.append(os.path.join(root, filename))

    #Read the lines for each file, if number of lines
    okTxtfiles=[]
    for txtfile in onlyFiles:
        with open(txtfile,'r') as f:
            if len(f.readlines())>0 or allAnnotations:
                okTxtfiles.append(os.path.split(txtfile)[-1])
    
    #remove start of string
    okTxtfiles=[txtfile[5:] if txtfile.startswith('data_') else txtfile for txtfile in okTxtfiles]
    
    #download file
    #Download in parallel
    executor = ThreadPoolExecutor(max_workers=10)
    futures = []    

    baseimagepath='http://roboweedsupport.com/Data' #'http://roboweedsupport.com/Data/233/20160518_135608.jpg'
    imageExt=['.jpg','.jpeg','.tif','.tiff','.png']
    #For all images, see if it exists, then download
    for txtfile in okTxtfiles:
        for extension in imageExt:
            

            if not os.path.exists(os.path.join(uploadid,txtfile.replace('.txt',extension))):
                #Copy from rws server, if computer is tbrain2
                if os.path.exists(os.path.join('/mnt/rwsdata/',uploadid,txtfile.replace('.txt',extension))):
                    #else download from server
                    a = executor.submit(os.system, 'cp "'+os.path.join('/mnt/rwsdata/',uploadid,txtfile.replace('.txt',extension))+'" "'+os.path.join(uploadid,txtfile.replace('.txt',extension))+'"')
                    futures.append(a)      
                else:
                    a = executor.submit(getImage, baseimagepath,uploadid,txtfile.replace('.txt',extension))
                    futures.append(a)            
    
    executor.shutdown(wait=True)    


#download an images from a website to a directory
def getImage(baseimagepath,Outimagepath,filename):
    teststr=baseimagepath+'/'+str(Outimagepath)+'/'+filename
    if not os.path.exists(str(Outimagepath)+'/'+filename): #if image has not allready been downloaded
        r = requests.head(teststr)
    
        if r.status_code==200: #200 = request accepted
            
            createIfNotExist(str(Outimagepath))
            cmd='wget -O "'+str(Outimagepath)+'/'+filename+'" "'+teststr+'"'
            os.system(cmd)

#If roboweedsupport has not been updated, use alternative download method
#This method simply looks through all possible imagenames for all upload ids. It is slow, but works
def getRoboWeedSupportUploadHack(uploadid):
    from BeautifulSoup import BeautifulSoup


    #Get a list of all images
    bs = BeautifulSoup(urllib2.urlopen("http://roboweedsupport.com/Administration/Images").read())
    table = bs.find(lambda tag: tag.name=='table' and tag.has_key('id') and tag['id']=="ctl00_ContentPlaceHolder1_TableImages") 
    rows = table.findAll(lambda tag: tag.name=='tr')
    imagelist=[]
    for row in rows:
        cols = row.findAll('td')
        cols = [ele.text.strip() for ele in cols]
        imagelist.append([ele for ele in cols if ele]) # Get rid of empty values
    
    #Only keep the image names
    imagelist=[x for x in imagelist if len(x)>3]
    imagelist=[x[2].encode('utf8') for x in imagelist]
    
    
    baseimagepath='http://roboweedsupport.com/Data' #'http://roboweedsupport.com/Data/233/20160518_135608.jpg'

    #Download in parallel
    executor = ThreadPoolExecutor(max_workers=100)
    futures = []    

    #For all images, see if it exists, then download
    for filename in imagelist:
        a = executor.submit(getImage, baseimagepath,uploadid,filename)
        futures.append(a)            
    
    executor.shutdown(wait=True)
    

def getAnnotationDictionary(annotations, speciesDictionary):
    #Remove empty lines
    annotations=[annotations[j].replace('\n','').replace('\r','') for j in range(len(annotations))]
    annotations=[annotations[j] for j in range(len(annotations)) if not annotations[j]==""]

    #get plant id and annotation
    plantidIdx=[j for j in range(len(annotations)) if annotations[j].find('PlantId: ')>=0]

    #Index where annotations starts
    annotationIdx=[j for j in range(len(annotations)) if annotations[j].find('[')>=0]
    
    #Index of ares
    areaIdx=[j for j in range(len(annotations)) if annotations[j].find('Area:')>=0]

    
    #Create annotation dict
    annotationDict=[]
    for ii in range(len(plantidIdx)):
        annotationDict_tmp={}
        #replacePlantid with species
        annotationDict_tmp['species']=speciesDictionary[''.join(re.findall("[+-]?\d+", annotations[plantidIdx[ii]]))]
        annotationDict_tmp['coordinates']=annotations[annotationIdx[ii]]
        if len(areaIdx)>0:
           # annotationDict_tmp['area']=annotations[areaIdx[ii]]
            annotationDict_tmp['area']=float(re.findall("[+-]?\d+\.?\d*", annotations[areaIdx[ii]])[0])
            
        annotationDict.append(annotationDict_tmp)
    return annotationDict


def exportThumbnails(annotationDict,dirname,imagepath):
    ExifThreadExecutor = ThreadPoolExecutor(max_workers=40)
    ExifFutures = []
    
    
    #load the image
    image = Image.open(imagepath)
    exifdata = {
        PIL.ExifTags.TAGS[k]: v
        for k, v in image._getexif().items()
        if k in PIL.ExifTags.TAGS
    }
     
    image=np.array(image)

    #Stretch HSV
    stretchHSV=False
    if stretchHSV:
        image=strechHSV(image)
        
    currentDir=os.path.dirname(os.path.realpath(__file__))
    
    for plant in annotationDict:
        #_, _, Ithumbnail =cropImageFromAnnotation(plant['coordinates'], image, annotateSquare=False)
        Ithumbnail, _, _ =cropImageFromAnnotation(plant['coordinates'], image, annotateSquare=True) 
        
        
        #Create directory if it does not exist
        iii=0
        
        #Increment iii until the name does not exist
        annotationfile,extension=os.path.splitext(os.path.basename(imagepath))
        #writefilename = os.path.join('export',dirname,plant['species'],str(iii)+'.tiff')
        def makeWriteFileName(iii):
            return os.path.join('export',dirname,plant['species'],annotationfile+'_'+str(iii)+'.tiff')
            #return os.path.join('export',dirname,plant['species'],str(iii)+'.tiff')
        
        writefilename = makeWriteFileName(iii)
        
        while os.path.exists(writefilename):
            iii+=1
            #writefilename = os.path.join('export',dirname,plant['species'],str(iii)+'.tiff')
            writefilename = makeWriteFileName(iii)



        #scipy.misc.imsave(os.path.join('export',dirname,annotations[plantidIdx[ii]],str(iii)+'.png'),Ialpha)





#EXPORT AS PNG
        ''' #Save as PNG
        meta = PngImagePlugin.PngInfo()
        
#                # copy metadata into new object
#                for k,v in PILimage.info.iteritems():
#                    meta.add_text(k, v, 0)


        #Add original exif data from input jpg
        for k,v in exifdata.iteritems():
            meta.add_text(k, str(v), 0)

        #add resolution (pixels per mm to the metadata)
        if 'area' in plant.keys():
            if float(plant['area'])<0: #if key is present, but area is not defined, asume frame (0.25m2) has been used
                plant['area']=0.5*0.5*(1.0*max(image.shape[:2])/min(image.shape[:2]))
            #convert image area to pixels per mm.
            pixelspermm = 1.0 / (np.sqrt(plant['area']*1000000 / (1.0*image.shape[0] * image.shape[1])))
            

            meta.add_text('pixelsPerMM',str(pixelspermm))
        else:
            print('EXPORT DATA IS OLD. CONSIDER REEXPORT FROM http://roboweedsupport.com/Administration/ViewUpload?uploadid='+uploadid)


        #Save the image
        PILimage.save(os.path.join('export',dirname,plant['species'],str(iii)+'.png'),pnginfo=meta)
        '''

        pixelspermm=None
        if 'area' in plant.keys():
            if float(plant['area'])<0.00025: #if area is not defined
                #See if area has been detected previously                       
                pixelsPermm = get_PixelsPermmtxt_indirOrParentDir(writefilename)
                if pixelsPermm:
                    plant['area']=image.shape[0]*image.shape[0]/(1000000*pixelsPermm*pixelsPermm)

                #If area still not found, look for barcode
                if float(plant['area'])<0.00025: #if area is not defined (less than 5x5cm)
                    #Try to look for barcode
                    print('FEJL: Nogle gange bliver areal INF. for eksempelv for upload 252')
                    print('file for QR:',dirname)
                    QRarea = QRreader()
                    QRarea.filepath=dirname #Use the directory, not single file name -> checks all images in folder
                    area, pixelsPrMeter = QRarea.getImageArea()
                    if pixelsPrMeter is not None:
                        pixelspermm=1.0*pixelsPrMeter/1000
                    
                    if pixelsPrMeter<1 or pixelsPrMeter>100000: #If area is unrealistic (remember i meters)
                        area=None
                    if area>0:
                        plant['area']=area
                    elif area is None: 
                        #if key is present, but area is still not defined, asume frame (0.25m2) has been used
                        plant['area']=0.5*0.5*(1.0*max(image.shape[:2])/min(image.shape[:2]))
                        
            #convert image area to pixels per mm.
            if pixelspermm is None:
                pixelspermm = 1.0 / (np.sqrt(plant['area']*1000000 / (1.0*image.shape[0] * image.shape[1])))                

        #Area is not defined (Only old uploads), asume 50x50cm frame
        if pixelspermm is None or 'area' not in plant.keys():
            plant['area']=0.5*0.5*(1.0*max(image.shape[:2])/min(image.shape[:2]))
            pixelspermm = 1.0 / (np.sqrt(plant['area']*1000000 / (1.0*image.shape[0] * image.shape[1])))    

        #Export txt-file with image resolution, if it does not exist
        if not os.path.exists(os.path.join('export',dirname,'pixelsPermm.txt')):
            exportPixelPermm(os.path.join('export',dirname),pixelspermm)

        #Test: save as tiff using adobe-deflate (lossless). Fylder mindre end png og kan have exifdata
        print('Exporting '+ writefilename)
        #Save the image (Use PIL, so that the quality can be set)
        PILimage=PIL.Image.fromarray(Ithumbnail)
        PILimage.info["Comment"] = imagepath
        createIfNotExist(os.path.join('export',dirname,plant['species']))
        PILimage.save(writefilename,compression="tiff_adobe_deflate")
        
        cmd1="exiftool -q -m -overwrite_original -TagsFromFile " + '"' + imagepath + '"' + ' '+ '"' + writefilename + '"'
        cmd2='exiftool -config "'+os.path.join(currentDir,'ExifTool_config').encode('utf8') + '" -q -m -overwrite_original -OriginalImage=' + '"' + imagepath + '"' + ' '+ '"' + writefilename+ '"'
        cmd3='exiftool -config "'+os.path.join(currentDir,'ExifTool_config').encode('utf8') + '" -q -m -overwrite_original -PixelsPerMm=' + '"' + str(pixelspermm) + '"' + ' '+ '"' + writefilename + '"'
        
        #added 2016-10-26
        cmd23='exiftool -config "'+os.path.join(currentDir,'ExifTool_config').encode('utf8') +'" -q -m -overwrite_original -OriginalImage=' + '"' + imagepath + '"' +  ' -PixelsPerMm=' + '"' + str(pixelspermm) + '"'  + ' '+ '"' + writefilename+ '"'

        
        #a = ExifThreadExecutor.submit(os.system, cmd1 + ' ; ' + cmd2 + ' ; ' + cmd3)
        a = ExifThreadExecutor.submit(os.system, cmd1 + ' ; ' + cmd23)        
        ExifFutures.append(a)            
    
    ExifThreadExecutor.shutdown(wait=True)


def strechHSV(image):
    print('HSV IS BEING STRECHED!')
    from skimage import color, exposure
#    from scipy import misc
    
    imageHSV=color.rgb2hsv(image)
    #p2, p98 = np.percentile(imageHSV[:,:,0], (1, 99))
    #imageHSV[:,:,0] = exposure.rescale_intensity(imageHSV[:,:,0], in_range=(p2, p98))
    p2, p98 = np.percentile(imageHSV[:,:,1], (1, 99))
    imageHSV[:,:,1] = exposure.rescale_intensity(imageHSV[:,:,1], in_range=(p2, p98))
    p2, p98 = np.percentile(imageHSV[:,:,2], (1, 99))
    imageHSV[:,:,2] = exposure.rescale_intensity(imageHSV[:,:,2], in_range=(p2, p98))
    image=(color.hsv2rgb(imageHSV)*255).astype(np.uint8)
    

    #p2, p98 = np.percentile(img, (2, 98))
    #image=exposure.rescale_intensity(image, in_range=(p2, p98))
    #misc.imsave('testimg.jpg',image)
    return image



def get_PixelsPermmtxt_indirOrParentDir(directory):
    temppath=directory
    #Look for pixelsPermm file in directory or parent directory, when found. read it and scale image acordingly
    while len(os.path.split(temppath)[0])>0:
        path=os.path.split(temppath)[0]
        if os.path.exists(os.path.join(path,'pixelsPermm.txt')):
            with open(os.path.join(path,'pixelsPermm.txt'),'r') as f:
                return float(f.readline())
            break
        else:
            temppath=os.path.split(temppath)[0]  




def exportPixelPermm(exportdir,pixelsPermm):   
        #if float(area)>0: #If area is defined
        createIfNotExist(exportdir)
        with open(os.path.join(exportdir,'pixelsPermm.txt'),'w') as f:
            f.write(str(pixelsPermm))

#Export the annotation data for a single object in image
#Denne funktion organiserer billederne, så de kan bruges til træning af RCNN
#Se https://github.com/zeyuanxy/fast-rcnn/blob/master/help/train/README.md 
#og
#http://sunshineatnoon.github.io/Train-fast-rcnn-model-on-imagenet-without-matlab/
def exportLocalizationDataforFastRCNN(imagepath,annotations,skipUnsorted=False):
    if skipUnsorted:
        ignoreLabels=['Unsorted','Junk']
        for label in ignoreLabels:
            annotations=[an for an in annotations if label not in an]
    if len(annotations)>0:
        #Run through all annotations
        labels=[]
        minXYmaxXY=[]
        
        jsonreader=json.JSONDecoder()
        
        for plant in annotations:
        
            jsonDict=jsonreader.decode(plant['coordinates'])
        
            #Get the x and y coordinates.
            xcoordinates=np.round(np.array([jsonDict[i]['X'] for i in range(len(jsonDict))])).astype(np.uint16)
            ycoordinates=np.round(np.array([jsonDict[i]['Y'] for i in range(len(jsonDict))])).astype(np.uint16)
            
            
            #Get the bounding box in the following format xmin, ymin, xmax, ymax
            xmin=min(xcoordinates)
            ymin=min(ycoordinates)        
            xmax=max(xcoordinates)
            ymax=max(ycoordinates)            
        
            #Get label
            labels.append(plant['species'].replace('æ','ae').replace('ø','oe').replace('å','aa').replace('Æ','Ae').replace('Ø','Oe').replace('Å','Aa'))
            minXYmaxXY.append([xmin,ymin,xmax,ymax])
    
        #get filename
        annotationfile=os.path.splitext(os.path.basename(imagepath))[0]
        annotationDir,filename=os.path.split(imagepath)
        
        #Save file with annotations
        AnnotationFullFileBegin="""<? w_819.xml ?>
<annotation>
  <folder>RWS</folder>
    <filename>
      <item>ITEM</item>
    </filename>
    """
    
        AnnotationFullFileEnd="""</annotation>"""

        Object="""    <object>
      <name>NAME</name>
      <bndbox>
      <xmin>XMIN</xmin>
      <ymin>YMIN</ymin>
      <xmax>XMAX</xmax>
      <ymax>YMAX</ymax>
      </bndbox>
    </object>
"""
        AnnotationStringList=[AnnotationFullFileBegin.replace('ITEM',filename)]
        #Create string to save
        for i in range(len(labels)):
            Objecttmp=Object
            Objecttmp=Objecttmp.replace('NAME',labels[i])
            Objecttmp=Objecttmp.replace('XMIN',str(minXYmaxXY[i][0]))
            Objecttmp=Objecttmp.replace('YMIN',str(minXYmaxXY[i][1]))
            Objecttmp=Objecttmp.replace('XMAX',str(minXYmaxXY[i][2]))
            Objecttmp=Objecttmp.replace('YMAX',str(minXYmaxXY[i][3]))
            AnnotationStringList.append(Objecttmp)
        
        AnnotationStringList.append(AnnotationFullFileEnd)
        if not os.path.exists(os.path.join('RCNNAnnotations',annotationDir,'Annotations')):
            os.makedirs(os.path.join('RCNNAnnotations',annotationDir,'Annotations'))
        with open(os.path.join('RCNNAnnotations',annotationDir,'Annotations',annotationfile+'.xml'), 'w') as f:
            f.writelines(AnnotationStringList)
  

#Export localization data for detectnet (Kitti format)                       
def exportDetectionDataKittiFormat(imagepath,annotations,skipUnsorted=False):
    '''
    #Values    Name      Description
    ----------------------------------------------------------------------------
       1    type         Describes the type of object: 'Car', 'Van', 'Truck',
                         'Pedestrian', 'Person_sitting', 'Cyclist', 'Tram',
                         'Misc' or 'DontCare'
       1    truncated    Float from 0 (non-truncated) to 1 (truncated), where
                         truncated refers to the object leaving image boundaries
       1    occluded     Integer (0,1,2,3) indicating occlusion state:
                         0 = fully visible, 1 = partly occluded
                         2 = largely occluded, 3 = unknown
       1    alpha        Observation angle of object, ranging [-pi..pi]
       4    bbox         2D bounding box of object in the image (0-based index):
                         contains left, top, right, bottom pixel coordinates
       3    dimensions   3D object dimensions: height, width, length (in meters)
       3    location     3D object location x,y,z in camera coordinates (in meters)
       1    rotation_y   Rotation ry around Y-axis in camera coordinates [-pi..pi]
       1    score        Only for results: Float, indicating confidence in
                         detection, needed for p/r curves, higher is better.    
                         
                         
      eg.
      Car    0.00  0   -2.12 749.24 184.57 946.06 282.98 0.002  0.02  0.01  3.87  1.64  12.49 -1.83
      type   trunc oc  alpha xmin   ymin   xmax   ymax   dim    dim   dim   loc   loc   loc   roty
    '''
    currentDir=os.path.dirname(os.path.realpath(__file__))
    image = np.array(Image.open(imagepath))
    imagesize=image.shape
    
    
    if skipUnsorted:
        ignoreLabels=['Unsorted','Junk']
        for label in ignoreLabels:
            annotations=[an for an in annotations if label not in an['species']]    
    if len(annotations)>0:
        
        #Run through all annotations        
        jsonreader=json.JSONDecoder()
        AnnotationStringList=[]
        for plant in annotations:
        
            jsonDict=jsonreader.decode(plant['coordinates'])
        
            #Get the x and y coordinates.
            xcoordinates=np.round(np.array([jsonDict[i]['X'] for i in range(len(jsonDict))])).astype(np.uint16)
            ycoordinates=np.round(np.array([jsonDict[i]['Y'] for i in range(len(jsonDict))])).astype(np.uint16)
    
            if 'BrushSize' in jsonDict[0]:
                if jsonDict[0]['BrushSize'] is not None:
                    brushsize=int(jsonDict[0]['BrushSize'])
                else:
                    brushsize=50 #Default value for old images, where brushsize was not yet implemented on IGIS' server
            else:
                brushsize=50 #Default value for old images, where brushsize was not yet implemented on IGIS' server
            
            #Get the bounding box in the following format xmin, ymin, xmax, ymax
            xmin=int(min(xcoordinates)-np.floor(brushsize/2))
            ymin=int(min(ycoordinates)-np.floor(brushsize/2))  
            xmax=int(max(xcoordinates)+np.floor(brushsize/2))
            ymax=int(max(ycoordinates)+np.floor(brushsize/2))  
            
                      
            
            #Get truncated area
            
            #Number of pixels of annotation that is in the image
            nPixelsInX=sum([1 for x in range(xmin,xmax) if 0<=x<=imagesize[1]])
            nPixelsInY=sum([1 for x in range(ymin,ymax) if 0<=x<=imagesize[0]])
            pixelsOfAnnotationInImage = nPixelsInX*nPixelsInY
            
            #size of boundingbox
            annotationarea = 1.0*(xmax-xmin)*(ymax-ymin)
            annotationOutOfImage = annotationarea - pixelsOfAnnotationInImage
            trunkation = annotationOutOfImage / annotationarea         
            
            trunkation=round(trunkation,2)
            


            #Kitti format
            if ' ' in plant['species']:
                plant['species']='"'+plant['species']+'"'
            AnnotationStringList.append([plant['species'],str(trunkation),'3','0',str(xmin),str(ymin),str(xmax),str(ymax),'0.03','0.0','0.0','0.0','0.0','0.0','0.0'])
            
            
    
        #get filename
        annotationfile,extension=os.path.splitext(os.path.basename(imagepath))
        annotationDir,filename=os.path.split(imagepath)
        
        #Save file with annotations


        if not os.path.exists(os.path.join('KittiAnnotations',annotationDir,'Annotations')):
            os.makedirs(os.path.join('KittiAnnotations',annotationDir,'Annotations'))
        with open(os.path.join('KittiAnnotations',annotationDir,'Annotations',annotationfile+'.txt'), 'w') as f:
            #f.writelines(AnnotationStringList)
            #f.writelines([' '.join(an)+'\n' for an in AnnotationStringList2])
            f.writelines([' '.join(an)+'\n' for an in AnnotationStringList])


        exportImageAndLabelsToSingleDirectory=True
        if exportImageAndLabelsToSingleDirectory:
            createIfNotExist(os.path.join('KittiAnnotationsOneDir','images'))  
            createIfNotExist(os.path.join('KittiAnnotationsOneDir','labels'))
            iiiii=0
            while os.path.exists(os.path.join('KittiAnnotationsOneDir','labels',str(iiiii)+'.txt')):
                iiiii+=1
            os.system('cp "'+imagepath+'" KittiAnnotationsOneDir/images/'+str(iiiii)+extension)
    
            with open(os.path.join('KittiAnnotationsOneDir','labels',str(iiiii)+'.txt'), 'w') as f:
                #f.writelines(AnnotationStringList)   
                f.writelines([' '.join(an)+'\n' for an in AnnotationStringList])


            createIfNotExist(os.path.join('KittiAnnotationsOneDir','labelsOneLabel'))
            with open(os.path.join('KittiAnnotationsOneDir','labelsOneLabel',str(iiiii)+'.txt'), 'w') as f:
                for i in range(len(AnnotationStringList)):
                    AnnotationStringList[i][0]='plant'
                #AnnotationStringList=[' '.join(an)+'\n' for an in AnnotationStringList]
                f.writelines([' '.join(an)+'\n' for an in AnnotationStringList])
           
            #Keep reference to original image
            cmd2='exiftool -config "'+os.path.join(currentDir,'ExifTool_config').encode('utf8') +'" -m -overwrite_original -OriginalImage=' + '"' + imagepath + '"' + ' '+ '"' + 'KittiAnnotationsOneDir/images/'+str(iiiii)+extension + '"'
            try:
                t2.join() #If thread is running, wait for it to finish before overwriting it
            except:
                pass
            t2 = Thread(target=os.system, args=(cmd2,))
            t2.start()


#Export localization data for YOLO                       
def exportLocalizationDataForYOLO(imagepath,annotations,skipUnsorted=False):
    if skipUnsorted:
        ignoreLabels=['Unsorted','Junk']
        for label in ignoreLabels:
            annotations=[an for an in annotations if label not in an['species']]  
    if len(annotations)>0:
            
        #Get original image size
        image = Image.open(imagepath)
        imsize=image.size
        del image
        
        #Run through all annotations
        labels=[]
        minXYmaxXY=[]
        
        jsonreader=json.JSONDecoder()
        
        for plant in annotations:
        
            jsonDict=jsonreader.decode(plant['coordinates'])
        
            #Get the x and y coordinates.
            xcoordinates=np.round(np.array([jsonDict[i]['X'] for i in range(len(jsonDict))])).astype(np.uint16)
            ycoordinates=np.round(np.array([jsonDict[i]['Y'] for i in range(len(jsonDict))])).astype(np.uint16)
    
            
            #Get the bounding box in the following format xmin, ymin, xmax, ymax
            xmin=min(xcoordinates)
            ymin=min(ycoordinates)        
            xmax=max(xcoordinates)
            ymax=max(ycoordinates)         
        
            #Convert to yolo standard, which is relative to image size
            X=((1.0*xmax-xmin)/2+xmin)/imsize[0]
            Y=((1.0*ymax-ymin)/2+ymin)/imsize[1]
            WIDTH=(1.0*xmax-xmin)/imsize[0]
            HEIGHT=(1.0*xmax-xmin)/imsize[1]
            
        
            #Get label
            #labels.append(plant['species'].replace('æ','ae').replace('ø','oe').replace('å','aa').replace('Æ','Ae').replace('Ø','Oe').replace('Å','Aa'))
            labels.append(plant['species'])
            minXYmaxXY.append([X,Y,WIDTH,HEIGHT])
    
        #get filename
        annotationfile=os.path.splitext(os.path.basename(imagepath))[0]
        annotationDir,filename=os.path.split(imagepath)
        
        #Save file with annotations
        
        Object='"NAME" X Y WIDTH HEIGHT\n' #x,y is center
    
        AnnotationStringList=[]
        #Create string to save
        for i in range(len(labels)):
            Objecttmp=Object
            Objecttmp=Objecttmp.replace('NAME',labels[i])
            Objecttmp=Objecttmp.replace('X',str(minXYmaxXY[i][0]))
            Objecttmp=Objecttmp.replace('Y',str(minXYmaxXY[i][1]))
            Objecttmp=Objecttmp.replace('WIDTH',str(minXYmaxXY[i][2]))
            Objecttmp=Objecttmp.replace('HEIGHT',str(minXYmaxXY[i][3]))
            AnnotationStringList.append(Objecttmp)
        
    
        if not os.path.exists(os.path.join('YOLOAnnotations',annotationDir,'Annotations')):
            os.makedirs(os.path.join('YOLOAnnotations',annotationDir,'Annotations'))
        with open(os.path.join('YOLOAnnotations',annotationDir,'Annotations',annotationfile+'.txt'), 'w') as f:
            f.writelines(AnnotationStringList)
            
  
            
            
def cropImageFromAnnotation(jsonentry, image, annotateSquare=False):
    jsonDict=jsonreader.decode(jsonentry)
    
    #Get the x and y coordinates.
    X=np.round(np.array([jsonDict[i]['X'] for i in range(len(jsonDict))])).astype(np.uint16)
    Y=np.round(np.array([jsonDict[i]['Y'] for i in range(len(jsonDict))])).astype(np.uint16)

    #brushsize=np.round(np.array([jsonDict[i]['brushsize'] for i in range(len(jsonDict))])).astype(np.uint16)
    if 'BrushSize' in jsonDict[0]:
        if jsonDict[0]['BrushSize'] is not None:              
            brushsize=int(jsonDict[0]['BrushSize'])
        else:
            print('Brushsize not defined, default to 50')
            brushsize=50 #Default value for old images, where brushsize was not yet implemented on IGIS' server
    else:
        brushsize=50 #Default value for old images, where brushsize was not yet implemented on IGIS' server
    brushsize=brushsize-np.mod(brushsize-1,2) #for even numbers, allways round to the odd number below


    #Create bounding box for plant
    Xmin=np.min(X)
    Xmax=np.max(X)
    Ymin=np.min(Y)
    Ymax=np.max(Y)


    cropSquared=True
    
    if cropSquared:
        #Find max of height and width
        sidelength=max([Xmax-Xmin,Ymax-Ymin])
        xcenter=(Xmax-Xmin)/2+Xmin
        ycenter=(Ymax-Ymin)/2+Ymin
        
        Xmin=np.floor(xcenter-sidelength/2)
        Ymin=np.floor(ycenter-sidelength/2)
        Xmax=np.ceil(xcenter+sidelength/2)
        Ymax=np.ceil(ycenter+sidelength/2)
    

    #If the brush is larger than the minimum index, the cropping will result in negative indices. 
    #Therefore, make sure that we are within the image
    Xmin=max(brushsize,Xmin)
    Ymin=max(brushsize,Ymin)
    Xmax=min(image.shape[1]-brushsize,Xmax)
    Ymax=min(image.shape[0]-brushsize,Ymax)

    #Make sure max is greater than min, if that is not the case, just make it a single pixel
    Xmax=max(Xmin,Xmax)
    Ymax=max(Ymin,Ymax)    
    
    Ymax=Ymax.astype(np.uint32)
    Xmax=Xmax.astype(np.uint32)
    Ymin=Ymin.astype(np.uint32)
    Xmin=Xmin.astype(np.uint32)
    
    Icropped = image[int((Ymin-np.floor(brushsize/2))-1):int(Ymax+np.floor(brushsize/2)),int(Xmin-np.floor(brushsize/2)-1):int(Xmax+np.floor(brushsize/2)),:]
  
    #If annotateSquare is true, we ignore the brush stroke and just color the whole alpha channel
    if not annotateSquare:
        #Create mask with same size as line. Brushsize is not considered yet
        try:
            Imask=np.zeros((Ymax-Ymin+1,Xmax-Xmin+1),dtype=np.uint8)
        except:
            pass
            
        #Set the annoations
        for i in range(len(X)):
            #Subtract one as Python index starts from 0
            try:
                Imask[Y[i]-Ymin,X[i]-Xmin]=1
            except:
                pass
    
        #Create structure alements from brush size
        strel=cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (brushsize,brushsize))
        
        #Allocate output from dilation. Make space for strel
        ImaskDilate=np.zeros(Imask.shape+np.floor(brushsize/2)*2,dtype=np.uint8)   
        
        #Put Imask in center of ImaskDilate
        ImaskDilate[np.ceil(brushsize/2):np.ceil(brushsize/2)+Imask.shape[0],np.ceil(brushsize/2):np.ceil(brushsize/2)+Imask.shape[1]]=Imask
        
        #Apply dilation
        ImaskDilate=cv2.dilate(ImaskDilate,strel)
    else:
        #ImaskDilate=np.ones(((Ymax-Ymin+1,Xmax-Xmin+1)+np.floor(brushsize/2)*2).astype(np.int),dtype=np.uint8)
        ImaskDilate=np.ones(Icropped.shape[:2],dtype=np.uint8)
        
    Imasked = Icropped*np.tile(ImaskDilate,(3,1,1)).transpose(1,2,0)
    Ialpha = np.concatenate((Icropped,np.expand_dims(ImaskDilate*255,2)),axis=2)

        
    
   # print(time.time()-t)
    return Icropped, Imasked, Ialpha
    


def unzipFile(url,path=''):
    zfile = zipfile.ZipFile(url)
    for name in zfile.namelist():
        (dirname, filename) = os.path.split(name)
        if not os.path.exists('./' + path+'/'+name):
            print("Decompressing " + filename + " on " + dirname)            
            #os.makedirs(dirname)
            zfile.extract(name,path)
    
def getZip(url):
    file_name = url.split('/')[-1]
    u = urllib2.urlopen(url)
    f = open(file_name, 'wb')
    meta = u.info()
    file_size = int(meta.getheaders("Content-Length")[0])
    print("Downloading: %s Bytes: %s" % (file_name, file_size))
    
    file_size_dl = 0
    block_sz = 8192
    while True:
        buffer = u.read(block_sz)
        if not buffer:
            break
    
        file_size_dl += len(buffer)
        f.write(buffer)
        status = r"%10d  [%3.2f%%]" % (file_size_dl, file_size_dl * 100. / file_size)
        status = status + chr(8)*(len(status)+1)
        #print status,
        #print status,"           \r",
        print(status, end='\r')    
    f.close()
        
if __name__ == '__main__':
    # Arguments - integers, except for binary/continous. Default uses binary.
    # Run with option --continuous for continuous output.
    import argparse
    parser = argparse.ArgumentParser(description='Command line options')    
    #parser.set_defaults(binary=True)    
    parser.add_argument('--uploadid', type=str, dest='uploadid', default='256')
    parser.add_argument('--skipUnsorted', action="store_true", dest='skipUnsorted')
    args = parser.parse_args(sys.argv[1:])
    main(**{k:v for (k,v) in vars(args).items() if v is not None})
