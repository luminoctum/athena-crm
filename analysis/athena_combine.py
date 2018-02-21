#! /usr/bin/env python
import argparse, os, shutil
from glob import glob
from datetime import datetime
from numpy import sort
from subprocess import check_call

## This function requies:
## 1) ncrcat
## 2) mppnccombine
## 3) cdo
def CombineNetcdfFiles(path, outfile = '', no_remove = False):
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
  for field in fields:
    for i in range(nb):
      print('Processing field "%s" block "%d" ...' % (field, i))
      files = path + '/' + '%s.block%d.%s.*.nc' % (case, i, field)
      target = '%s.%s.nc.%04d' % (outfile, field, i)
      check_call('ncrcat -h %s -o %s' % (files, target), shell = True)
      check_call('ncatted -O -a %s,%s,%c,%c,%d %s' %
        ('NumFilesInSet', 'global', 'c', 'i', nb, target), shell = True)

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
    check_call('cdo merge %s %s.nc' % (singles, outfile), shell = True)
    for f in singles.split():
      os.remove(f)
  else:
    shutil.move(outfile, '%s.%s.nc' % (outfile, field), '%s.nc' % outfile)
  print('Finish case "%s", output file is %s.nc' % (case, outfile))

if __name__ == '__main__':
  now = datetime.now()
  today = '%02d%02d' % (now.month, now.day),
  files = sort(glob('log.%s?' % today))
  if len(files) > 1:
    latest = files[-1].split('.')[-1]
  else:
    latest = today + 'a'
    check_call('touch log.%s' % latest)

  parser = argparse.ArgumentParser()
  parser.add_argument('-d', '--dir',
    default = '.',
    help = 'directory of the simulation to combine'
    )
  parser.add_argument('-o', '--output',
    default = latest,
    help = 'appending additional name to the output'
    )
  parser.add_argument('-n', '--no-remove',
    action = 'store_true',
    help = 'do not remove original files'
    )
  args = vars(parser.parse_args())

  CombineNetcdfFiles(args['dir'], outfile = args['output'], no_remove = args['no_remove'])
