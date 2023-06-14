"""所有要動的參數都在這裡"""
LA= 5                     # 指施力點的高程(3m,5m,7m)
I = 595210.94	 * 0.01**4  # 慣性矩
A = 1387.09	 / 10000      # 截面積

#直接勁度法
from numpy import matrix, linalg, sqrt , power, cos, sin, pi , arctan
P = 180000  #kgf
L=7 #m
LB=7-LA #m
E = 2.04*10**10  #kg/cm2 * (cm/m)^2 = kg/m2
RA = P*(LB**2)*(3*LA+LB)/(L**3)  #kgf
RB = P*(LA**2)*(3*LB+LA)/(L**3)  #kgf
MA = P*LA*LB**2/L**2  #kgf*m
MB = P*LB*LA**2/L**2  #kgf*m
La = 7
Lb = sqrt(53)

ka_local = matrix([
    [E*A/La,0,0,-E*A/La,0,0],
    [0,12*E*I/La**3,-6*E*I/La**2,0,-12*E*I/La**3,-6*E*I/La**2],
    [0,-6*E*I/La**2,4*E*I/La,0,6*E*I/La**2,2*E*I/La],
    [-E*A/La,0,0,E*A/La,0,0],
    [0,-12*E*I/La**3,6*E*I/La**2,0,12*E*I/La**3,6*E*I/La**2],
    [0,-6*E*I/La**2,2*E*I/La,0,6*E*I/La**2,4*E*I/La],
])

kb_local = matrix([
    [E*A/Lb,0,0,-E*A/Lb,0,0],
    [0,12*E*I/Lb**3,-6*E*I/Lb**2,0,-12*E*I/Lb**3,-6*E*I/Lb**2],
    [0,-6*E*I/Lb**2,4*E*I/Lb,0,6*E*I/Lb**2,2*E*I/Lb],
    [-E*A/Lb,0,0,E*A/Lb,0,0],
    [0,-12*E*I/Lb**3,6*E*I/Lb**2,0,12*E*I/Lb**3,6*E*I/Lb**2],
    [0,-6*E*I/Lb**2,2*E*I/Lb,0,6*E*I/Lb**2,4*E*I/Lb],
])

angle_a =  pi/2
angle_b =  -arctan (7/2)
Te_a = matrix([
    [cos(angle_a),sin(angle_a),0,0,0,0],
    [-sin(angle_a),cos(angle_a),0,0,0,0],
    [0,0,1,0,0,0],
    [0,0,0,cos(angle_a),sin(angle_a),0],
    [0,0,0,-sin(angle_a),cos(angle_a),0],
    [0,0,0,0,0,1]
])

Te_b =matrix([
    [cos(angle_b),sin(angle_b),0,0,0,0],
    [-sin(angle_b),cos(angle_b),0,0,0,0],
    [0,0,1,0,0,0],
    [0,0,0,cos(angle_b),sin(angle_b),0],
    [0,0,0,-sin(angle_b),cos(angle_b),0],
    [0,0,0,0,0,1]
])

ka_globe = matrix.transpose(Te_a) * ka_local * Te_a
kb_globe = matrix.transpose(Te_b) * kb_local * Te_b

Γa = matrix([
    [1,0,0,0,0,0,0,0,0],
    [0,1,0,0,0,0,0,0,0],
    [0,0,1,0,0,0,0,0,0],
    [0,0,0,1,0,0,0,0,0],
    [0,0,0,0,1,0,0,0,0],
    [0,0,0,0,0,1,0,0,0]
])

Γb = matrix([
    [0,0,0,1,0,0,0,0,0],
    [0,0,0,0,1,0,0,0,0],
    [0,0,0,0,0,1,0,0,0],
    [0,0,0,0,0,0,1,0,0],
    [0,0,0,0,0,0,0,1,0],
    [0,0,0,0,0,0,0,0,1]
])

Ka = matrix.transpose(Γa) * ka_globe * Γa
Kb = matrix.transpose(Γb) * kb_globe * Γb
K_system = Ka + Kb

K_need = K_system[3:6,3:6]
R_test = matrix([
    [RB],
    [0],
    [-MB]
])
r_test = linalg.inv(K_need) * R_test

from numpy import concatenate
zero1 = matrix([
    [0],
    [0],
    [0]
])

zero2 = matrix([
    [0],
    [0],
    [0]
])

r_test_full = concatenate((zero1,r_test,zero2), axis = 0)

FEM_test = matrix([
    [0],
    [RA],
    [-MA],
    [0],
    [RB],
    [MB]
])

Qa = ka_local * Te_a * Γa *r_test_full  + FEM_test
Qb = kb_local * Te_b * Γb *r_test_full
print(f'''頂點水平位移為{r_test[0]}m
頂點垂直位移為{r_test[1]}m
頂點轉角為    {r_test[2]}rad
''')
print(f'''垂直桿件的內力
r1={Qa[0]}kgf
r2={Qa[1]}kgf
r3={Qa[2]}kgf*m
r4={Qa[3]}kgf
r5={Qa[4]}kgf
r6={Qa[5]}kgf*m
      ''')
print(f'''斜撐桿件的內力
r1={Qb[0]}kgf
r2={Qb[1]}kgf
r3={Qb[2]}kgf*m
r4={Qb[3]}kgf
r5={Qb[4]}kgf
r6={Qb[5]}kgf*m
      ''')

#能量法計算位移
G = 0.78 *10**10

#垂直桿件
nLN_EA1 = Qa[0]**2 * L / (180000 * E * A )
vV_GA1 = Qa[1]**2 * LA / (180000*G*A) +  Qa[4]**2 * LA / (180000*G*A)
mM_EI1 = (  LA / (180000 * E * I) * (Qa[2]**2 + Qa[2]*Qa[1]*LA + Qa[1]**2 * LA **2 / 3 ) +
            LB / (180000 * E * I) * (Qa[5]**2 - Qa[5]*Qa[4]*LA + Qa[4]**2 * LA **2 / 3 ) )

#斜桿件
nLN_EA2 = Qb[0]**2 * Lb / (180000 * E * A )
vV_GA2 = Qa[1]**2 * Lb / (180000*G*A)
mM_EI2 = Lb / (180000 * E * I) * (Qb[2]**2 + Qb[2]*Qb[1]*Lb + Qb[1]**2 * Lb **2 / 3 )

print("施力點變位為",nLN_EA1 + nLN_EA2 + vV_GA1 + vV_GA2 + mM_EI1 + mM_EI2 ,'m')

Fy = 2.530            # tf/cm2
E  = 2.04 * 10**3     # tf/cm2
L1 = 7 *100           # cm
L2 = sqrt(53) * 100   # cm
k  = 1
A  = A * 10000              # cm2
#Pu = abs(Qa[0])/1000            # tf
#M2 = abs(Qa[2]/10)      # tf * cm
#M1 = abs(Qa[5]/10)       # tf * cm

rx = 20.8            # cm
ry = 11.9              # cm
Zx = 27.2   *10**3       # cm3
Sx = 21   *10**3
bf = 455 /10          # cm
d  = 569 /10
tf = 125 /10
rt = 144 /10            # cm

#ASD法

Af = tf * bf

#算軸力強度
Cc = sqrt(2*pi**2*E/Fy)
print("Cc=",Cc)

kL1r = max (k*L1/rx , k*L1/ry)
kL2r = max (k*L2/rx , k*L2/ry)

Lc = min(20 * bf / sqrt(Fy), 1400/(d*Fy/(Af)))

'''----------------------檢驗直桿--------------------------'''

#軸壓強度

if kL1r > Cc :
    print("彈性挫曲")

else:
    print("非彈性挫曲")
    R1 = kL1r / Cc
    Fa1 = ( 1 - R1**2/2) / (5/3 + 3/8 *R1 -1*R1**3 /8) * Fy
    print('Fa1 = ',Fa1)

#算彎矩強度

if abs(Qa[2]) > abs(Qa[5]):
  M2 = Qa[2]/10
  M1 = Qa[5]/10
  print('M2=',M2,'M1=',M1)
else:
  M1 = Qa[2]/10
  M2 = Qa[5]/10
  print('M2=',M2,'M1=',M1)

Cb1 = 1   #Cb,Cm不共存

Lu1 = rt*sqrt(7170*Cb1/Fy)

Le1 = rt*sqrt(35800*Cb1/Fy)

print("Lc=",Lc,"cm","Lu1=",Lu1,"Le1=",Le1)

if L1 < Lc:
  Fb1 = 0.66 * Fy

elif Lc < L1 and L1 < Lu1 :
  Fb1 = 0.6 * Fy

else:
  Fb1 = max(2/3 - Fy*(L1/rt)**2/(107600*Cb1), 840*Cb1/(L1*d/Af))*Fy
  if Fb1 > 0.6*Fy:
    Fb1 = 0.6*Fy

print('Fb1=',Fb1)


#直桿檢驗
fb1_top = abs(Qa[2]/10) / Sx
fb1_bottom = abs(Qa[5]/10) / Sx
fa1 = abs(Qa[0]/1000) / A
print(f'''fa1 ='{fa1}
fb1_top = {fb1_top}
fb1_bottom = {fb1_bottom}
      ''')

Cm = 0.6 - 0.4*(M1/M2)     #雙曲率
Fex = (12/23)* pi**2 * E / (L1/rx)**2

if fa1/Fa1 > 0.15:
  print(f'''
直桿
  大軸力分析：
  頂部
  穩定性準則={(fa1) / Fa1 + Cm * fb1_top / ( (1-fa1/Fex)*Fb1 )}
  降伏性準則={fa1/(0.6*Fy) + fb1_top / Fb1}
  底部
  穩定性準則={(fa1) / Fa1 + Cm * fb1_bottom / ( (1-fa1/Fex)*Fb1 )}
  降伏性準則={fa1/(0.6*Fy) + fb1_bottom / Fb1}

''')

if fa1/Fa1 < 0.15:
  print(f'''
直桿
  頂部
  小軸力準則={fa1/(Fa1) + fb1_top / Fb1}
  底部
  小軸力準則={fa1/(Fa1) + fb1_bottom / Fb1}
''')


'''----------------------檢驗斜桿--------------------------'''

#斜桿軸壓強度

if kL2r > Cc :
    print("彈性挫曲")

else:
    print("非彈性挫曲")
    R2 = kL2r / Cc
    Fa2 = ( 1 - R2**2/2) / (5/3 + 3/8 *R2 -1*R2**3 /8) * Fy
    print('Fa2 = ',Fa2)

#斜桿彎矩強度

if abs(Qb[2]) > abs(Qb[5]):
  M2 = Qb[2]/10
  M1 = Qb[5]/10
  print('M2=',M2,'M1=',M1)
else:
  M1 = Qb[2]/10
  M2 = Qb[5]/10
  print('M2=',M2,'M1=',M1)

Cb2 = 1   #Cb,Cm不共存

Lu2 = rt*sqrt(7170*Cb2/Fy)

Le2 = rt*sqrt(35800*Cb2/Fy)

print("Lc=",Lc,"cm","Lu=",Lu2,"Le=",Le2)

if L2 < Lc:
  Fb2 = 0.66 * Fy

elif Lc < L2 and L2 < Lu2 :
  Fb2 = 0.6 * Fy

else:
  Fb2 = max(2/3 - Fy*(L2/rt)**2/(107600*Cb2), 840*Cb2/(L2*d/Af))*Fy
  if Fb2 > 0.6*Fy:
    Fb2 = 0.6*Fy

print('Fb2=',Fb2)

#斜桿檢驗
fb2_top = abs(Qb[2]/10) / Sx
fb2_bottom = abs(Qb[5]/10) / Sx
fa2 = abs(Qb[0]/1000) / A

Cm = 0.6 - 0.4*(M1/M2)     #雙曲率
Fex = (12/23)* pi**2 * E / (L2/rx)**2

if fa2/Fa2 > 0.15:
  print(f'''
斜桿
  大軸力分析：
  頂部
  穩定性準則={(fa2) / Fa2 + Cm * fb2_top / ( (1-fa2/Fex)*Fb2 )}
  降伏性準則={fa2/(0.6*Fy) + fb2_top / Fb2}
  底部
  穩定性準則={(fa2) / Fa2 + Cm * fb2_bottom / ( (1-fa2/Fex)*Fb2 )}
  降伏性準則={fa2/(0.6*Fy) + fb2_bottom / Fb2}

''')

if fa2/Fa2 < 0.15:
  print(f'''
斜桿
  頂部
  小軸力準則={fa2/(Fa2) + fb2_top / Fb2}
  底部
  小軸力準則={fa2/(Fa2) + fb2_bottom / Fb2}
''')

