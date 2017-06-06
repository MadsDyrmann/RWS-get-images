# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 11:50:13 2016

@author: mads

#Make crops of images and annotations


History:
    2017-05-23: Added support for keeping  both original and cropped
    2017-05-24: Fixed an error causing multiple threads to access same image when writing
"""

import os
import csv
import codecs
from PIL import Image
import numpy as np
from concurrent.futures import ThreadPoolExecutor
executor = ThreadPoolExecutor(max_workers=30)
exportExecutor = ThreadPoolExecutor(max_workers=30)

futures = []
import time



def createIfNotExist(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)    
        
        
def flipAndRotateImageAndAnnotations(image,annotations,interval):
    #interval is number between 0 and 7
    if interval<4:
        image = image.transpose(Image.TRANSPOSE)
        for i in range(len(annotations)):
            annotations[i][4:8]=[annotations[i][ii] for ii in [5,4,7,6]]
            
    for i in range(np.mod(interval,4)):
        annotations = rotateCoordinate90degrees(annotations,imagesize=image.size)
        image = image.transpose(Image.ROTATE_90)

    return image, annotations


def rotateCoordinate90degrees(annotations,imagesize):
    annotationsold=np.copy(annotations).tolist()
    for ix,annotation in enumerate(annotations):
        annotations[ix][5]=str(int(imagesize[0])-int(annotationsold[ix][4]))
        annotations[ix][4]=annotationsold[ix][5]
        annotations[ix][7]=str(int(imagesize[0])-int(annotationsold[ix][6]))
        annotations[ix][6]=annotationsold[ix][7]

    return annotations



def getAnnotationsInBbox(annotationList, ystart,xstart,(cropsizeX,cropsizeY)):
    #Run through all annotations and convert them
    AnnotationStringList=[]

    for annotation in annotationList:
        label=annotation[0]
        xmin=int(float(annotation[4]))-xstart
        ymin=int(float(annotation[5]))-ystart
        xmax=int(float(annotation[6]))-xstart
        ymax=int(float(annotation[7]))-ystart
        
        
        replacePlantBycar=False
        if replacePlantBycar:
            if label.lower()=='plant':
                label='Car'
        
        #Get truncated area
        
        #Number of pixels of annotation that is in the image
        nPixelsInX=sum([1 for x in range(xmin,xmax) if 0<=x<=cropsizeX])
        nPixelsInY=sum([1 for x in range(ymin,ymax) if 0<=x<=cropsizeY])
        pixelsOfAnnotationInImage = nPixelsInX*nPixelsInY
        
        #size of boundingbox
        annotationarea = 1.0*(xmax-xmin)*(ymax-ymin)
        annotationOutOfImage = annotationarea - pixelsOfAnnotationInImage
        trunkation = annotationOutOfImage / annotationarea
                   
        
        #truncate crop
        xmin=max(0,xmin)
        ymin=max(0,ymin)
        xmax=min(cropsizeX,xmax)
        ymax=min(cropsizeY,ymax)
        
        
        #Get the size of the annotations after cropping
        annotationarea = 1.0*(xmax-xmin)*(ymax-ymin)
        

        # object size should be between 50×50 and 400×400 px in the input images
        #only accept annotations where at least 20% is visible (1.0-0.2 = 0.8)
        okAnnotations = trunkation<0.8
        #Only accept objects larger than 50 pixels on the widthest point
#            okAnnotations = okAnnotations and 50<(xmax-xmin) and 50<(ymax-ymin)
        #Only accept objects smaller than 400px on the widthest point
#            okAnnotations = okAnnotations and (xmax-xmin)<400 and (ymax-ymin)<400
        
        if okAnnotations:      #only accept annotations where at least 80% is visible (1.0-0.8 = 0.2)
            Object=annotation[:]
       
            #Create string to save
            Object[0]=label
            
            Object[1]=str(trunkation)
            Object[4]=str(xmin)
            Object[5]=str(ymin) #ændret til xmin ymin xmax ymax format
            Object[6]=str(xmax)
            Object[7]=str(ymax)
            AnnotationStringList.append(Object[:])
    return AnnotationStringList

def exportImageAndAnnotation(imcrop,AnnotationStringList, writeName='', flipAndRotate=False):
    images = []
    AnnotationStringListRot=[]   
       
    if flipAndRotate:
        intervals=[1,3,4,6] # number from 0 to 7
        for interval in intervals:
            imrot,anno = flipAndRotateImageAndAnnotations(imcrop.copy(),np.copy(AnnotationStringList).tolist(),interval=interval)
            images.append(imrot)
            AnnotationStringListRot.append(np.copy(anno).tolist())
    else:
        
        images.append(imcrop)
        AnnotationStringListRot=[AnnotationStringList]

    #Export the annotation
    createIfNotExist(os.path.join(annotationOutputPath,trainval))
    createIfNotExist(os.path.join(imageOutputPath,trainval))        

    iii=0

    for imcrop, AnnotationStringList in zip(images,AnnotationStringListRot):
        while os.path.exists(os.path.join(annotationOutputPath,trainval,writeName+str(iii)+'.txt')):
            iii+=1
        
        #MAXANNOTATIONS=-1
        MAXANNOTATIONS=50
        AnnotationStringList=AnnotationStringList[:MAXANNOTATIONS]
        
        
        annotationsavename = os.path.join(annotationOutputPath,trainval,writeName+str(iii)+'.txt')
        imagesavename = os.path.join(imageOutputPath,trainval,writeName + str(iii)+image_extension)

        with open(annotationsavename, 'w') as f:
            f.writelines([' '.join(an)+'\n' for an in AnnotationStringList])
            
        imcrop.save(imagesavename)
        
        assert os.path.exists(imagesavename) and os.path.exists(annotationsavename)

        
        #Copy exifdata
        cmd1="exiftool -overwrite_original -TagsFromFile " + '"' + imagefile + '"' + ' '+ '"' + imagesavename + '"'
  
        a = executor.submit(os.system, cmd1)
        futures.append(a)




###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################










imagedir='/media/mads/Eksternt drev/GetThumbnailsFromServer/KittiAnnotationsOneDir/images'
annotationPath='/media/mads/Eksternt drev/GetThumbnailsFromServer/KittiAnnotationsOneDir/labelsOneLabel'

imageOutputPath='/media/mads/Eksternt drev/GetThumbnailsFromServer/KittiAnnotationsOneDir/ImagesCrop'
annotationOutputPath='/media/mads/Eksternt drev/GetThumbnailsFromServer/KittiAnnotationsOneDir/LabelsCrop'


#imagedir='/media/slow-storage/DataFolder/MadsDyrmann/Database/GetThumbnailsFromServerForDetection/KittiAnnotationsOneDir/images'
#annotationPath='/media/slow-storage/DataFolder/MadsDyrmann/Database/GetThumbnailsFromServerForDetection/KittiAnnotationsOneDir/labelsOneLabel'
#
#imageOutputPath='/media/slow-storage/DataFolder/MadsDyrmann/Database/GetThumbnailsFromServerForDetection/KittiAnnotationsOneDir/ImagesCrop'
#annotationOutputPath='/media/slow-storage/DataFolder/MadsDyrmann/Database/GetThumbnailsFromServerForDetection/KittiAnnotationsOneDir/LabelsCrop'












cropsize={}
#Kitti size
#cropsize['width']=1224
#cropsize['height']=370

cropsize['width']=1224
cropsize['height']=1024

#Probabilities for val/train/test (Nb: Digits does not support both test and val for object detection -> set test to 0.0)
Trainprob=0.9 #Probability of ending in train
TestProb=0.0 #val prob is the remainer
futures=[]


datestr=time.strftime("%Y-%m-%d")
imageOutputPath=os.path.join(imageOutputPath,'imagesCropped'+'_'+datestr+'_'+str(cropsize['height'])+'x'+str(cropsize['width']))
annotationOutputPath=os.path.join(annotationOutputPath,'labelsDontCareCropped'+'_'+datestr+'_'+str(cropsize['height'])+'x'+str(cropsize['width']))



imagefiles = []
for root, dirnames, filenames in os.walk(imagedir):
    for filename in filenames:
        if filename.lower().endswith(('.png','.jpg','.jpeg','.tif','.tiff','.gif')):
            imagefiles.append(os.path.join(root, filename))
imagefiles.sort()    


#imagefiles=imagefiles[:5]



iii=0
for i in range(len(imagefiles)):
    print('image '+str(i)+' of ' + str(len(imagefiles))+', '+imagefiles[i])
    imagefile=imagefiles[i]
    #get image name
    _,filename=os.path.split(imagefile)
    filename, image_extension = os.path.splitext(filename)
    image_extension=image_extension.lower()

    im = Image.open(imagefile)

    #Pad image if it is too small
    padw=max(cropsize['width']-im.size[0],0)
    padh=max(cropsize['height']-im.size[1],0)
    if padw>0 or padh>0:
        im=Image.fromarray(np.pad(im,((0,padh),(0,padw),(0,0)),'constant'))

    possibleStartHeight=im.height-cropsize['height']
    possibleStartWidth=im.width-cropsize['width']

    #ncrops=5
    ncrops=int(np.ceil(1.0*(im.size[0]*im.size[1])/(cropsize['width']*cropsize['height'])))
    
    rowStart=np.random.randint(0,possibleStartHeight+1,ncrops)
    colStart=np.random.randint(0,possibleStartWidth+1,ncrops)
    
    
    
    
    

    #Load the annotation
    speciesDictionary=[]
    with codecs.open(annotationPath+'/'+filename+'.txt','r') as csvfile:
        reader = csv.reader(csvfile,delimiter=' ',skipinitialspace=True)
        annotationList=[row for row in reader]

    #Use image for train, validation or test
    p=np.random.rand()
    
    if p<Trainprob:
        trainval='train'
    elif p>Trainprob and p<Trainprob+TestProb:
        trainval='test'
    else:
        trainval='val'
    
    
    #Run through all the crops
    for ystart,xstart in zip(rowStart,colStart):
                
        AnnotationStringList = getAnnotationsInBbox(annotationList, ystart,xstart,(cropsize['width'],cropsize['height']))
        
        
        if len(AnnotationStringList)>0:

            #Create cropped image
            imcrop = im.crop((xstart, ystart, xstart+cropsize['width'], ystart+cropsize['height']))
            imname=im.filename.split('/')[-1].split('.')[-2]
            
            #exportExecutor.submit(exportImageAndAnnotation,imcrop,AnnotationStringList, writeName=imname+'_'+str(ystart)+'_'+str(xstart),flipAndRotate=True)
            exportImageAndAnnotation(imcrop,AnnotationStringList, writeName=imname+'_'+str(ystart)+'_'+str(xstart)+'_',flipAndRotate=True)
            
            
    #Export whole image, scaled
    imcrop = im.resize((cropsize['width'],cropsize['height']))
    scalefactor = (1.0*cropsize['width']/im.size[0],1.0*cropsize['height']/im.size[1])
    AnnotationStringList = np.copy(annotationList).tolist()
    for i in range(len(AnnotationStringList)):
        AnnotationStringList[i][4] = str(int(int(AnnotationStringList[i][4])*scalefactor[0]))
        AnnotationStringList[i][5] = str(int(int(AnnotationStringList[i][5])*scalefactor[1]))
        AnnotationStringList[i][6] = str(int(int(AnnotationStringList[i][6])*scalefactor[0]))
        AnnotationStringList[i][7] = str(int(int(AnnotationStringList[i][7])*scalefactor[1]))
    
    
    #exportExecutor.submit(exportImageAndAnnotation,imcrop,AnnotationStringList, writeName=imname,flipAndRotate=True)
    exportImageAndAnnotation(imcrop,AnnotationStringList,flipAndRotate=True)

    while exportExecutor._work_queue.full(): #Avoid filling the ram, by not reading more images into memory than what fits in the work queue
        pass
    
executor.shutdown(wait=True)
exportExecutor.shutdown(wait=True)












#import matplotlib.pyplot as plt
#import seaborn as sns
#sns.set(style="darkgrid")
#
## Load the long-form example gammas dataset
#gammas = sns.load_dataset("gammas")
#
## Plot the response with standard error
#plt.imshow(x)
#sns.tsplot(data=gammas, time="timepoint", unit="subject", condition="ROI", value="BOLD signal")
#
#
#plt.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='off')
#plt.tick_params(axis='both', which='both', bottom='off', top='off', labelbottom='off')
#plt.legend().set_visible(False)
#plt.show()