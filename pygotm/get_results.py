from pygotm import medsea

def get_var(varname='temp',hour_range=range(0,365*24),year=2013,run='ASM3-75m'):
     from numpy.ma import masked_all
     from time import time
     fn = 'results-{:s}-{:d}.nc'.format(run,year)
     yrdays = 366 if year%4 == 0 else 365
     dim = (yrdays*24,122,medsea.M,medsea.N) if varname == 'temp' or varname == 'salt' else (yrdays*24,medsea.M,medsea.N)
     var = masked_all(dim)     
     elapsed = 0
     count = 0
     start = time()
     for m,n in medsea.sea_mn:
         local_folder = medsea.get_local_folder(m,n)
         join = medsea.os.path.join
         fullfn = join(local_folder,fn)
         isfile = medsea.os.path.isfile
         if not(isfile(fullfn)):
             continue
         with medsea.Dataset(fullfn,'r') as ds:
              if ds[varname].shape[0] != dim[0]:
                   continue
              if varname == 'temp' or varname == 'salt':
                   var[:,:,m,n] = ds[varname][hour_range,::-1,0,0]
              else:
                   var[:,m,n] = ds[varname][hour_range,0,0]
         elapsed = time()-start
         count += 1
         if count%100 == 0:
             C = medsea.sea_m.size
             est = elapsed/count*(C-count)
             print('{:d} / {:d} grid points, {:.1f} s elapsed, {:02d} m {:.1f} s to go.'.format(count,C,elapsed,int(est/60),est%60), end='\r')
     return var
