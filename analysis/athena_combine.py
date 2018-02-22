#! /usr/bin/env python
import argparse, os, shutil
from glob import glob
from datetime import datetime
from numpy import sort
from subprocess import check_call
from multiprocessing import Pool

def CombineSingleBlock(arglist):
  path, case, field, outfile, nb, i = arglist
  print('Processing field "%s" block "%d" ...' % (field, i))
  files = path + '/' + '%s.block%d.%s.*.nc' % (case, i, field)
  target = '%s.%s.nc.%04d' % (outfile, field, i)
  check_call('ncrcat -h %s -o %s' % (files, target), shell = True)
  check_call('ncatted -O -a %s,%s,%c,%c,%d %s' %
      ('NumFilesInSet', 'global', 'c', 'i', nb, target), shell = True)

## This function requies:
## 1) ncrcat, ncks
## 2) mppnccombine
def CombineNetcdfFiles(path, outfile = '', no_remove = False, threads = 1):
  files = glob(path + '/*.out*.[0-9][0-9][0-9][0-9][0-9].nc')
  cases = []
  blocks = []
  fields = []
  stamps = []

  for fname in files:
    case, block, field, stamp, ext = os.path.basename(fname).split('.')
    if case not in cases:
      cases.append(case)
    if block not in blocks:
      blocks.append(block)
    if field not in fields:
      fields.append(field)
    if stamp not in stamps:
      stamps.append(stamp)

  if len(cases) > 1:
    assert False, 'More than one case to combine'
  else:
    case = cases[0]

  nb, nf, ns = len(blocks), len(fields), len(stamps)
  if len(files) != nb*nf*ns:
    assert False, 'Missing fiels'
  else:
    print('Number of blocks: %d' % nb)
    print('Number of fields: %d' % nf)
    print('Number of stamps: %d' % ns)

  if outfile == '':
    outfile = case
  else:
    outfile = case + '-' + outfile

  singles = ''
  fields = sorted(fields, reverse = True)
  pool = Pool(threads)
  for field in fields:
    arglist = []
    for i in range(nb):
      arglist.append((path, case, field, outfile, nb, i))
    pool.map(CombineSingleBlock, arglist)

    print('Combining blocks ...')
    if not os.path.exists('mppnccombine'):
      check_call('gcc -O3 -o mppnccombine ../analysis/mppnccombine.c -lnetcdf',
        shell = True)
    check_call('./mppnccombine %s.%s.nc' % (outfile, field), shell = True)

    print('Removing temporary files ...')
    for f in glob('%s.%s.nc.????' % (outfile, field)):
      os.remove(f)

    if not no_remove:
      print('Removing blocks ...')
      for f in glob(path + '/' + '%s.block*.%s.*.nc' % (case, field)):
        os.remove(f)
      print('Done for field "%s".' % field)

    singles += '%s.%s.nc ' % (outfile, field)

  if len(fields) > 1:
    print('Combining %s' % singles)
    check_call('ncks -A %s' % singles, shell = True)
    shutil.move(singles.split()[-1], '%s.nc' % outfile)
  else:
    shutil.move('%s.%s.nc' % (outfile, field), '%s.nc' % outfile)

  for f in singles.split()[:-1]:
    os.remove(f)
  print('File processed: %s'% singles)
  print('Finish case "%s", output file is %s.nc' % (case, outfile))

def CleanNetcdfFiles():
  files1 = glob('*.out?.nc.????')
  files2 = []
  for f in files1:
    print 'cleaning file:', f
    os.remove(f)
    if f[:-5] not in files2:
      files2.append(f[:-5])
  for f in files2:
    print 'cleaning file:', f
    os.remove(f)
  print 'Done.'

if __name__ == '__main__':
  now = datetime.now()
  today = '%02d%02d' % (now.month, now.day)

  parser = argparse.ArgumentParser()
  parser.add_argument('-d', '--dir',
    default = '.',
    help = 'directory of the simulation to combine'
    )
  parser.add_argument('-o', '--output',
    default = today,
    help = 'appending additional name to the output'
    )
  parser.add_argument('-n', '--no-remove',
    action = 'store_true',
    help = 'do not remove original files'
    )
  parser.add_argument('-c', '--clean',
    action = 'store_true',
    help = 'clean temporary files'
    )
  parser.add_argument('-p', '--threads',
    default = '4',
    help = 'number of parallel threads to use'
    )
  args = vars(parser.parse_args())

  if args['clean']:
    CleanNetcdfFiles()
  else:
    CombineNetcdfFiles(args['dir'], outfile = args['output'], 
        no_remove = args['no_remove'], threads = int(args['threads']))
