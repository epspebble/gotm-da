def get_4D_var(varname='temp',hour_range=range(0,31*24),outfn='results-ASM3-75m-2013.nc'):
     from numpy.ma import masked_all
     from time import time
     var = masked_all((len(hour_range),122,medsea.M,medsea.N))
     elapsed = 0
     count = 0
     start = time()
     for m,n in medsea.sea_mn:
         local_folder = medsea.get_local_folder(m,n)
         join = medsea.os.path.join
         fn = join(local_folder,outfn)
         isfile = medsea.os.path.isfile
         if not(isfile(fn)):
             continue
         with medsea.Dataset(fn,'r') as ds:
             var[:,:,m,n] = ds[varname][hour_range,::-1,0,0]
         elapsed = time()-start
         count += 1
         if count%100 == 0:
             C = medsea.sea_m.size
             est = elapsed/count*(C-count)
             print('{:d} / {:d} grid points, {:.1f} s elapsed, {:02d} m {:.1f} s to go.'.format(count,C,elapsed,int(est/60),est%60), end='\r')
     return var

def get_4D_var(varname='temp',hour_range=range(0,31*24),outfn='results-ASM3-75m-2013.nc'):
     from numpy.ma import masked_all
     from time import time
     var = masked_all((len(hour_range),122,medsea.M,medsea.N))
     elapsed = 0
     count = 0
     start = time()
     for m,n in medsea.sea_mn:
         local_folder = medsea.get_local_folder(m,n)
         join = medsea.os.path.join
         fn = join(local_folder,outfn)
         isfile = medsea.os.path.isfile
         if not(isfile(fn)):
             continue
         with medsea.Dataset(fn,'r') as ds:
             var[:,:,m,n] = ds[varname][hour_range,::-1,0,0]
         elapsed = time()-start
         count += 1
         if count%100 == 0:
             C = medsea.sea_m.size
             est = elapsed/count*(C-count)
             print('{:d} / {:d} grid points, {:.1f} s elapsed, {:02d} m {:.1f} s to go.'.format(count,C,elapsed,int(est/60),est%60), end='\r')
     return var

def get_3D_var(varname='x-taus',outfn='results-ASM3-75m-2013.nc'):
     from numpy.ma import masked_all
     from time import time
     var = masked_all((8760,medsea.M,medsea.N))
     elapsed = 0
     count = 0
     medsea.tic()
     start = time()
     for m,n in medsea.sea_mn:
         local_folder = medsea.get_local_folder(m,n)
         join = medsea.os.path.join
         fn = join(local_folder,outfn)
         isfile = medsea.os.path.isfile
         if not(isfile(fn)):
             continue
         with medsea.Dataset(fn,'r') as ds:
             var[:,m,n] = ds[varname][:,0,0]
         elapsed = time()-start
         count += 1
         if count%100 == 0:
             C = medsea.sea_m.size
             est = elapsed/count*(C-count)
             print('{:d} / {:d} grid points, {:.1f} s elapsed, {:02d} m {:.1f} s to go.'.format(count,C,elapsed,int(est/60),est%60), end='\r')
     return var