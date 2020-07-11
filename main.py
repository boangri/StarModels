# import pandas as pd
import phys as ph
from Sun import SSM18

print('speed of light = %e' % ph.c)
print('Solar mass =  %e' % SSM18.M)

df = SSM18.load_interpolated_data()

print(df)
print(ph.Pressure(100, 13e7, 0.732, 1, 0.02))