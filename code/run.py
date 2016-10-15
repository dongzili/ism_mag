import spec as sp

data = sp.read()
#sp.plot_spec(data,save='y')

ic = sp.sec_spec(data)
sp.plot_spec(ic,save='y',spec_type='sec')
