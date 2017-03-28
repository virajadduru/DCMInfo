
# -*- coding: utf-8 -*-
"""
This is my module for extracting dicom information from datasets.

@author: vradduru
"""
#%% imports modules

import dicom as dcmtool
import os
import cPickle
import subprocess

import numpy as np
import logging as log
import sys
import shutil
import time

#%% FUNCTION DEFINITIONS

def get_dcm_params(dcm_files,dcm_field_names,data_dir=''):
    """
    Extracts the dicom fields from the dicom files supplied. It forces on 
    files that dont have dcm signature
    
    inputs:
    
        dcm_files: List of the dicom file names you want to process
        
        dcm_field_names: List of the fields you want to extract from the files
        
        data_dir: (default = '') source directory in case you want to specify
    
    output: 
        
        dcm_field_vals: A dict object containing dicom fields as keys and their 
            values in a list corresponding to the keys in the order of dcm_files.
    
    """
    num_files = len(dcm_files)
    dcm_field_vals=dict([[name,[]] for name in dcm_field_names])
    if type(dcm_files) == list:        
        for i in range(num_files):
            # read each file
            dcm_file = dcmtool.read_file(os.path.join(data_dir,dcm_files[i]),force = True)
            
            # append series description value of each file to series_desc
            for field in dcm_field_names:
                dcm_field_vals.get(field).append(dcm_file.get(field))
            
    else:
        assert type(dcm_files) == str
        dcm_file = dcmtool.read_file(os.path.join(data_dir,dcm_files),force = True)
        for field in dcm_field_names:
            dcm_field_vals.update({field:dcm_file.get(field)})
                
    dcm_field_vals.update({'FileNames':dcm_files})    
    return dcm_field_vals
    
    
    
def group_image_indices(field_vals_list):
    """
    Identifies the images in dicoms and gives out the list of indices belonging
    to each image
    
    Iteratively groups the indices starting with the primary field and subgroups the 
    primary groups using secondary field and viceversa
    
    Inputs: 
    
    Outputs:
    """ 
    
    field_sizes = [len(field_vals_list[k]) for k in range(len(field_vals_list))]
    
    field_sizes = list(set(field_sizes))
    
    assert len(field_sizes)==1
    
    field_size = field_sizes[0]
#    field_vals_list = [primary_field_vals]+list(secondary_field_vals)
    grouped_indices=[range(field_size)]
    
    for field_vals in field_vals_list:    
        
        assert len(field_vals)==field_size
        
        temp_grouped_indices = []
        for idx in grouped_indices:
            
            field_unique=[]    
            field_grouped_indices = [];
            
            for j in idx:
                try:
                    group_id = field_unique.index(field_vals[j])
                    field_grouped_indices[group_id].append(j)
                except:
                    field_unique.append(field_vals[j])
                    field_grouped_indices.append([j])
            
            temp_grouped_indices+=field_grouped_indices
        grouped_indices = temp_grouped_indices
            
    return grouped_indices
    
    
    
def group_field_vals(field_vals,grouped_indices):
    """
    field_vals must be a single list of values to be grouped
    grouped_indices must be a list of lists
    """
    grouped_field_vals = [[field_vals[grouped_indices[p][k]] \
                            for k in range(len(grouped_indices[p]))] \
                            for p in range(len(grouped_indices))] 
    
    return grouped_field_vals
   
   
def write_tocsv(data_list,filename):
    with open(filename,'wt') as csvfile:
        for i in data_list:
            str_list = [str(k) for k in i]
            csvfile.write(','.join(str_list)+'\n')



def copy_subdata(source_file,dest_file,copylog_filename = r'copy_log.txt'):
    dest_dir = os.path.split(dest_file)[0]
    copylog_file=file(os.path.join(dest_dir,copylog_filename),'at')
    command = ["copy",source_file,dest_file]
    try:
        print '...copying'
        print source_file,'to'
        print dest_file
#        print subprocess.check_output(command)
        shutil.copyfile(source_file,dest_file)
        copylog_file.write(str(command)+'\n')
    
    finally:        
        copylog_file.close()


def open_dcm_file(file_name,data_dir=''):
    file_path = os.path.join(data_dir,file_name)
    command = ["C:\Program Files\RadiAntViewer64bit\RadiAntViewer.exe",file_path]
    subprocess.Popen(command)

# both the functions below are for saving dicom info of slices of an image
def save_info(dcm_object_list,file_name):
    f = file(file_name,'wb')
    for obj in dcm_object_list:
        cPickle.dump(obj,f)
    f.close()

    
def sortnsave(dcm_object_list, sort_with, file_name = None):
    tosort = np.asarray(sort_with,dtype='float32')
    sortidx = np.argsort(tosort)
    sorted_obj = []
    for obj in dcm_object_list:
        sorted_obj.append(np.asarray(obj,dtype='float32')[sortidx])
    
    if file_name != None:
        save_info(sorted_obj,file_name)
    
    return sorted_obj,sortidx

    
def getNifti(dicom_dir,convert=False,sub_nifti_path=None,log_file=None,compress = True):
    """
    Saves compressed images
    """
    if compress == True:
        c_arg = '-z'
        nii_ext = '.nii.gz'
    else:
        c_arg = ''
        nii_ext = '.nii'
    
    if sub_nifti_path == None:
        sub_nifti_path = dicom_dir

    if os.path.isdir(sub_nifti_path) :
        temp_file_list = os.listdir(sub_nifti_path)
        nifti_files = [[k,os.path.getsize(os.path.join(sub_nifti_path,k))/1024.0] for k in temp_file_list if k.endswith(nii_ext)]
    else :
        nifti_files = []
        
    if len(nifti_files)==0 and convert==True:
        print('\n...converting to nifti')
        if not os.path.isdir(sub_nifti_path):
            os.makedirs(sub_nifti_path)
        try:
            conv_log_txt = subprocess.check_output([r"C:\Users\vradduru\Tools\mricrogl_20160930\dcm2niix.exe",'-m',c_arg,'-f','%s_%p','-o',sub_nifti_path,dicom_dir],
                                                   cwd = r'C:\Users\vradduru\Tools\mricrogl_20160930')
            print conv_log_txt
        except:
            #e = sys.exc_info()[0]
            er = sys.exc_info()[1]
            # display error on console
            print er
            conv_log_txt = 'Error: With nifti conversion'
            
        
        if log_file != None:
            try:    
                nifticonvlog_file=file(log_file,'at')
                nifticonvlog_file.write('\n\n************************************************************\n')
                nifticonvlog_file.write('subject: '+os.path.split(sub_nifti_path)[-1]+'\n')
                nifticonvlog_file.write(conv_log_txt+'\n')
                
            except:
                print 'Error: Failed to log nifti convertion'
            
            finally:
                nifticonvlog_file.close()
            
        temp_file_list = os.listdir(sub_nifti_path)
        nifti_files = [[k,os.path.getsize(os.path.join(sub_nifti_path,k))/1024.0] for k in temp_file_list if k.endswith(nii_ext)]
        
    return nifti_files

    
def dicomFolderInfo(dataDir,dcmFieldNames,groupFieldOrder = None, single_file = False):
    
    # get dcm file list
    files_list = os.listdir(dataDir)
    dcm_files = [k for k in files_list if k.endswith('.dcm')]
    
    # get dcm_field values for all files
    if single_file == False :
        dcm_field_vals = get_dcm_params(dcm_files,dcmFieldNames,dataDir)
        
        # grouping data
        if groupFieldOrder != None:          
            data_for_grouping = [dcm_field_vals.get(groupFieldOrder[k]) for k in range(len(groupFieldOrder))] # pack grouping fields for grouping
            series_indices = group_image_indices(data_for_grouping) # get indices in groups
            grouped_field_vals=dict([[name,group_field_vals(dcm_field_vals.get(name),series_indices)] for name in dcm_field_vals.keys()]) # group all the field values
    
        else:
            grouped_field_vals = dcm_field_vals
            
        
        return grouped_field_vals    
        
    else:
        dcm_field_vals = get_dcm_params(dcm_files[0],dcmFieldNames,dataDir)
        
        return dcm_field_vals
    
def listDicomFolders(datadir, verbose = False,checkdcmfiles = False):
    dcm_folders = []
    dir_gen = os.walk(datadir)
    count = 0
    for r in dir_gen :
    # Get list of .dcm files
        dcm_list = [k for k in r[2] if k.endswith('.dcm')]
        if len(dcm_list) != 0:
            count += 1
            if verbose:
                print count,r[0]
            dcm_folders.append([r[0],dcm_list])
    return dcm_folders 
    
    
def savesublist(data_dir,filename,verbose = False):
    
    rawdata_dir = data_dir
    
    dcm_folders = listDicomFolders(rawdata_dir, verbose)
    
    subjects_filename = filename
    
    dcm_subsession = []
    for i in dcm_folders:
        path,session = os.path.split(i[0])
        sub = os.path.split(path)[-1]
        dcm_subsession.append([sub,session])
    
    dcm_subsession.sort(key=lambda x: x[0])
    with open(subjects_filename,'wt') as subjects_file:
        for sub,session in dcm_subsession:
            name = sub+'_'+session
            print name
            subjects_file.write(name+'\n')
            
    

    
def myfunc(sub_name,**extra):
    
    subj_dir = sub_name.replace('_','\\')    
    data_dir = os.path.join(extra.get('dcm_dir'),subj_dir)
    
    imageinfo_filename=r'imageinfo.txt'
    
    manual_copy = extra.get('manual_copy')
    dest_dir = [extra.get('dest_dir')]
    copy_suffix = '.nii.gz'
    
    convert_to_nii = extra.get('convert_to_nii')#True
    nifti_dir = extra.get('nifti_dir')#r'G:\DATA_NIFTI\MR-Head-DeID\data'
    sub_nifti_dir = os.path.join(nifti_dir,sub_name)    
    conv_info_file = os.path.join(sub_nifti_dir,'log_nifticonv.txt')
    
    write_dcm_info = True
    forcewrite_dcminfo = extra.get('forcewrite_dcminfo')#False
    if forcewrite_dcminfo == None:
        forcewrite_dcminfo = False
    dcminfo_filename = r'dcminfo.txt'
    
    
    

    updateDCMInfo = extra.get('updateDCMInfo')
    

    
    print 'Subject: ',sub_name  
    dcminfo_filepath = os.path.join(nifti_dir,sub_name,dcminfo_filename)
    
    if os.path.isfile(dcminfo_filepath) and forcewrite_dcminfo == False:
        with file(dcminfo_filepath,'rt') as dcminfofile:
            dcminfo = eval(dcminfofile.readline())
        
        
    else:
        print '\n...reading dicoms'
        dcminfo = {'subjectName':sub_name}
        dcm_field_names = ['SeriesDescription','SeriesNumber','ImageType',
                           'SliceThickness','SliceLocation','SpacingBetweenSlices','PixelSpacing',
                           'InstanceNumber','AcquisitionNumber','ImagesInAcquisition',
                           'ImagePositionPatient','ImageOrientationPatient','GantryDetectorTilt']
        
        grouping_fieldorder = ['SeriesDescription','SeriesNumber','ImageType']
           
        tic = time.clock()
        grouped_field_vals = dicomFolderInfo(data_dir,dcm_field_names,grouping_fieldorder)
        toc = time.clock()
        
        series_indices = grouped_field_vals.get('FileNames')    
        tot_files = sum([len(k) for k in series_indices])
        
        print tot_files,'files in',toc-tic,'seconds'
    # temporary
        y = dcmtool.read_file(os.path.join(data_dir,series_indices[1][0]))
        dcminfo.update({'Manufacturer':str(y.Manufacturer)})
        dcminfo.update({'ManufacturerModelName':str(y.ManufacturerModelName)})
        dcminfo.update({'PatientSex':str(y.get('PatientSex'))})
        dcminfo.update({'PatientAge':str(y.get('PatientAge'))})
        #dcminfo.update({'MagneticFieldStrength':str(y.get('MagneticFieldStrength'))})
        

    
        dcminfo.update({'PatientBirthDate':str(y.get('PatientBirthDate'))})    
        dcminfo.update({'StudyDate':str(y.get('StudyDate'))}) 
        
        
    
        param_bunch = [list(set(grouped_field_vals.get('SeriesNumber')[k]))+ \
                        list(set(grouped_field_vals.get('SeriesDescription')[k]))+ \
                        [len(grouped_field_vals.get('SeriesNumber')[k])]+ \
                        list(set(grouped_field_vals.get('GantryDetectorTilt')[k])) + \
                        #list(set(grouped_field_vals.get('SliceThickness')[k])) + \
                        [grouped_field_vals.get('PixelSpacing')[k][0]] + \
                        [list(set(grouped_field_vals.get('SliceThickness')[k]))]# + \
                        #list(set(grouped_field_vals.get('SpacingBetweenSlices')[k])) # + \
                       # list(set(grouped_field_vals.get('SeriesDescription')[k])) \
                        for k in range(len(series_indices)) ]
                            
        dcminfo.update({'images':param_bunch})
        
        
    if updateDCMInfo!=None :
        print '\n...Updating dcminfo'
        ukeys = updateDCMInfo.keys()
        if len(ukeys)!=0:
            for ukey in ukeys:
                dcminfo.update({ukey : updateDCMInfo.get(ukey)})
                print 'DCM_Info updated with: ',ukey
        print '\n'
    else:
        ukeys = []

    print '\nPatient age: ',dcminfo.get('PatientAge')
    print 'Manufacturer:',dcminfo.get('Manufacturer')
    print 'Model Name:',dcminfo.get('ManufacturerModelName')
    print 'Dicom Series Info:'
    param_bunch = dcminfo.get('images')
    for m in range(len(param_bunch)): print param_bunch[m]
    
## get nifti information
    nifti_files = getNifti(dicom_dir = data_dir, convert = convert_to_nii, sub_nifti_path = sub_nifti_dir, log_file = conv_info_file)        
    num_nifti_files = len(nifti_files)   
    print '\n',num_nifti_files,'nifti files found:'
    print nifti_files
    
# writing dicom info to file.
    if  (write_dcm_info == True and os.path.isfile(dcminfo_filepath) == False) or len(ukeys)!=0 or forcewrite_dcminfo == True :
        print '\nDCM info:'
        print dcminfo        
        with file(dcminfo_filepath,'wt') as dcminfofile:
            dcminfofile.write(str(dcminfo))
        print '\nDicom info written to file: ',dcminfo_filepath    
            
            
    ## manual copy
    
    if manual_copy == True:
        print '\nCopying Nifti file'
        for k in range(num_nifti_files):
            print k+1,nifti_files[k]
        copy_idx = input('Enter index of file: ')-1
        
        
        if copy_idx in range(num_nifti_files):
            copyniftifile = nifti_files[copy_idx][0]
            copy_subdata(os.path.join(nifti_dir,sub_name,copyniftifile),
                       os.path.join(dest_dir[0],sub_name+copy_suffix))#nifti_files[copy_idx][0]))
            print '\nCopying image info'  
            
            file_name = os.path.splitext(os.path.splitext(copyniftifile)[0]  )[0]
            
            series_num = int(file_name.split('_')[0])
            seriesid_list = [int(k[0]) for k in param_bunch if k[0]!=None]
            
            param_id = [k for k in range(len(seriesid_list)) if seriesid_list[k]==series_num]

            
            if len(param_id)==0:
                print 'No matching image info found'
                for m in range(len(param_bunch)): print m+1,param_bunch[m]
                info_idx = input('Enter index of series: ')-1
                
            elif len(param_id)==1:
                info_idx = param_id[0]
                print '1 matching image info found'
                
            else:
                print len(param_id),'matching Image info found'
                for m in range(len(param_id)): print m+1,param_bunch[param_id[m]]
                info_idx = param_id[input('Enter index of series: ')-1]
                
                    
            try:
                imageinfo_filepath = os.path.join(dest_dir[0],imageinfo_filename)
                file_exists = os.path.isfile(imageinfo_filepath)
                
                imageinfo_file=file(imageinfo_filepath,'at')
                if not file_exists:
                    header = [['Patient_ID','PatientAge','PatientSex','SeriesNumber','SeriesDescription','NumberofSlices','GantryDetectorTilt','PixelSpacing','SliceThickness','Manufacturer','ManufacturerModelName']]
                    imageinfo_file.write(str(header)+'\n')
                imageinfo = str([sub_name]+[dcminfo.get('PatientAge')]+[dcminfo.get('PatientSex')]+param_bunch[info_idx]+[dcminfo.get('Manufacturer'),dcminfo.get('ManufacturerModelName')])
                print imageinfo
                imageinfo_file.write(imageinfo+'\n')
            
            finally:
                imageinfo_file.close()
        else:
            print 'Copying nothing from this subject!'




