# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import os

#%%

def segHist(img,seg_msk = None,plt_func=None,nbins = 500, hist_range=None,cmap = None, normed = False, disp = True):
    """
    Plots a histogram of the image intensities with different colors for each
    segmentation.
    
    Input:

        img: Image for which histogram should be plotted.
        
        seg_msk: (default: class 1 on all voxels) Segmentation mask containing classifications for each voxel.     
        
        plt_func: In case you want to specify a matplotlib handle.
    
        nbins: (default: 500) Custom number of bins can be specified
            
        hist_range: (default: min to max of image data) User can select custom intensity range for histogram  
        
        cmap: (default: matplotlib default) User can select a custom colormap  
        
        normed: normalizes frequency by TOTAL size of the image. # should we do this with respect to seg_mask??
    
    Output:
    
        Returns frequency and binedges for each class in a list.
        
    """
    import numpy as np
    import matplotlib.cm as cm
    
    img_data = img.flatten()
    if seg_msk is not None:
        seg_data = seg_msk.flatten()
    else :
        seg_data = np.ones(img_data.shape)
    
    temp = np.where(seg_data>0)
    img_data = img_data[temp]
    seg_data = seg_data[temp]
    
    classes = np.unique(seg_data)
    classes = classes[np.where(classes != 0)[0]]
    
    numClasses = len(classes)
    if numClasses < 3:
        numClasses = 3
    colormap = cm.get_cmap(name = cmap)
    classcolor= colormap(np.arange(1/(numClasses*2.),1,1/float(numClasses))) #['blue','green','red']
    
    if plt_func==None:
        import matplotlib.pyplot as plt_func
        
    if hist_range==None:
        hist_range = (img_data.min(),img_data.max())
        
    plt_func.hold(True)
    output = []
    for i in classes:
        n, binedgs = np.histogram(img_data[np.where(seg_data==i)],nbins,hist_range)
        if normed:
            n = n/float(img_data.size)
        plt_func.plot(binedgs[1:],n,color=classcolor[i-1])
        output.append((n,binedgs))
    plt_func.hold(False)   

    return output
    
    
    
def openNifti(data_dir,subfolders = False,filename = None,fileext = '.nii'):
    
    """
    Sequentially open nifti images in a datafolder for viewing in ITKSnap 
    """    
    
    import subprocess
    import logging as log
    import sys
    import os
    
    dir_content = os.listdir(data_dir)
    dir_content = [os.path.join(data_dir,k) for k in dir_content]
    
    if subfolders == True:
        assert filename != None
        subjdirs = [i for i in dir_content if os.path.isdir(i)]
        files_list = [os.path.join(i,filename) for i in subjdirs]
        nifti_files = [[k,os.path.getsize(k)/1024.0] for k in files_list if os.path.isfile(k)]
        
    else:
        nifti_files = [[k,os.path.getsize(os.path.join(data_dir,k))/1024.0] for k in dir_content if k.endswith(fileext)]
    
    nifti_files = [k for k in nifti_files if os.path.isfile(k[0])]
    
    if nifti_files:
        numFiles = len(nifti_files)
        print 'Found',numFiles,'files\n'
        count = 0;
        for i in nifti_files:
            count = count+1
            print count ,r'/', numFiles
            command = [r"C:\Program Files\ITK-SNAP 3.6\bin\ITK-SNAP.exe",i[0]]
            try:
                print 'Opening : ',i
                print subprocess.check_output(command)
            
            except:
                er = sys.exc_info()[1]
                log.error(er)
    
    else:
        print 'No files found'
    
    
    
    
#%%
    
if __name__=='__main__':
    
    import nibabel.nifti1 as niftitool
    import numpy as np

    with open(os.path.join(data_dir,'subjectlist.txt')) as f:
        subject_list = f.read().splitlines()
        

    numplots = 50
    want_subplots = False
    
    if want_subplots is True:
        plotnum = numplots
    else:
        plotnum = 1
        
    f,pltarr = plt.subplots(plotnum,sharex=True)
    
    plt.suptitle('Intensity histogram of T2 image across 51 miccai08 subjects')
    for n,i in enumerate(subject_list):
        print i
        mode_image_file = os.path.join(data_dir,i,
                                       i+'_T2_masked_roi_restore_std.nii.gz')
        seg_mask_file = os.path.join(data_dir2,i,
                                     'segmented',
                                     i+'_T1_brain_roi_pveseg.nii.gz')

        mode_image = niftitool.load(mode_image_file)
        seg_mask = niftitool.load(seg_mask_file)
        
        if isinstance(pltarr,np.ndarray):
            plotter = pltarr[n]
            plt.setp(plotter.get_yticklabels(),visible=False)
        else:
            plotter = pltarr
            
        plt.setp(plotter.get_yticklabels(),visible=False)
        segHist(mode_image.get_data(),seg_mask.get_data(),plotter)
        if n==numplots-1:
            break
        
    plt.xlabel('Tissue Intensity')