#!/usr/bin/env python-sirius

from tkinter import Tk
from tkinter import filedialog
from shutil import copyfile
import os
import sys

readme_template = [
'Sirius BC Dipole',
'================',
'',
'current official model version : MODEL_VERSION',
'',
'fieldmap analysis folder       : FIELDMAP_FOLDER',
'',
'fieldmap file                  : FIELDMAP_FILE',
''
'',
'',
'Observations',
'============',
'',
'- This folder contains copies of relevant files from the official magnet model',
'- each time a new magnet model version becomes the official version "update-official-model.py" has to be edited and run',
]

def select_rawfield_fieldmap_analysis():

    # we don't want a full GUI, so keep the root window from appearing
    Tk().withdraw()
    # show an "Open" dialog box and return the path to the selected file
    rawfield_filename = filedialog.askopenfilename(initialdir = "../",
                                         title = "Select fieldmap analysis rawfield input file",
                                         filetypes = (("text data files","rawfield.in"),("all files","*.*")))

    # --- rawfield_filename ---
    #rawfield_filename = '/home/fac_files/data/sirius/lnls-ima/bo-dipoles/model-09/analysis/fieldmap/analysis01/rawfield.in'
    if not rawfield_filename: return None

    # --- fmap_analysis_folder ---
    if '/' in rawfield_filename:
        words = rawfield_filename.split('/')
    else:
        words = rawfield_filename.split('\\')
    fmap_analysis_folder = words[-2]

    # --- magnet_model ---
    magnet_model_version = None
    for i in range(len(words)):
        if 'analysis' in words[i]:
            magnet_model_version = words[i-1]
            break
    if not magnet_model_version:
        for word in words:
            if 'model' in word:
                magnet_model_version = word
                break

    # --- fieldmap_filename ---
    with open(rawfield_filename) as f:
        content = f.readlines()
    for line in content:
        words = line.split()
        if words and words[0] == 'fmap_filename':
            (dirName, fieldmap_filename) = os.path.split(words[1].replace("'",''))

    r = {
        'rawfield_filename':rawfield_filename,
        'fmap_analysis_folder':fmap_analysis_folder,
        'magnet_model_version':magnet_model_version,
        'fieldmap_filename':fieldmap_filename}

    return r


def write_readme(r):

    def _template_subs(words_prev, words_next, text):
        for i in range(len(text)):
            for j in range(len(words_prev)):
                text[i] = text[i].replace(words_prev[j], words_next[j])
        return text

    _subs = {
        'MODEL_VERSION' : r['magnet_model_version'],
        'FIELDMAP_FOLDER' : r['fmap_analysis_folder'],
        'FIELDMAP_FILE' : r['fieldmap_filename'],
    }

    readme_txt = _template_subs(list(_subs.keys()), list(_subs.values()), readme_template)
    fp = open('README.md','w')
    for line in readme_txt:
        fp.write(line+'\n');



def copy_files(r):
    os.system('cp -raf ../' + r['magnet_model_version'] + '/documentation ./')
    os.system('cp -raf ../' + r['magnet_model_version'] + '/simulation/magnetic/' + r['fieldmap_filename'] + ' ./fieldmap-file.txt')
    os.system('cp -raf ../' + r['magnet_model_version'] + '/analysis/fieldmap/' + r['fmap_analysis_folder'] + '/analysis.txt ./')
    os.system('cp -raf ../' + r['magnet_model_version'] + '/analysis/fieldmap/' + r['fmap_analysis_folder'] + '/analysis.pdf ./')
    os.system('cp -raf ../' + r['magnet_model_version'] + '/analysis/fieldmap/' + r['fmap_analysis_folder'] + '/field_on_trajectory.txt ./')
    os.system('cp -raf ../' + r['magnet_model_version'] + '/analysis/fieldmap/' + r['fmap_analysis_folder'] + '/multipoles.txt ./')
    os.system('cp -raf ../' + r['magnet_model_version'] + '/analysis/fieldmap/' + r['fmap_analysis_folder'] + '/trajectory.txt ./')
    os.system('cp -raf ../' + r['magnet_model_version'] + '/analysis/fieldmap/' + r['fmap_analysis_folder'] + '/*.svg ./')

def main(argv):
    print('selecting analysed rawfield file...')
    r = select_rawfield_fieldmap_analysis()
    print(' rawfield_filename   : ', r['rawfield_filename'])
    print(' fmap_model_version  : ', r['fmap_analysis_folder'])
    print(' magnet_model_version: ', r['magnet_model_version'])
    print(' fieldmap_filename   : ', r['fieldmap_filename'])
    print('writting README.md file in this folder')
    write_readme(r)
    print('copying analysis files to this official folder')
    copy_files(r)

if __name__ == '__main__':
    main(sys.argv)
