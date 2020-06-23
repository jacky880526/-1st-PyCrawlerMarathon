import numpy as np 

s0 = 50
r = 0.1
q = 0.05
sigma = 0.8
T_t=0.25
t = 0.25
sim = 10000
rep = 20
n = 100
s_ave_t = 50
m= 100
K= 50

def Euro_call(s0, r, q, sigma, T_t, t, n, s_ave_t, m, K):
	dt = (T_t)/n
	u = np.exp(sigma * np.sqrt(dt))
	d = 1/u
	p = (np.exp((r-q)*(dt))-d) / (u-d)	
	b_n = (n/T_t) * t      #n before t
	b_s = s_ave_t*(b_n+1)  
	
	st = np.zeros((n+1,n+1))
	s_ave = np.zeros((n+1,n+1,m+1))
	call_value = np.zeros((n+1,n+1,m+1))
	
	st[0,0] = s0
	for j in range(0,n+1):           #n-period
		st[n,j] = s0*(u**(n-j))*(d**j)
		if n-j == j:
			st[n,j] = s0
	
	for j in range(0,n):
		st[n-1,j] = s0*(u**(n-1-j))*(d**j)
		if n-1-j == j:
			ST[n-1,j] = s0

	for i in range(1,n-1):
		for j in range(0,i+1): 
			if (n-i) % 2 == 0:
				st[i,j] = st[n,j+int((n-i)/2)]
			else:
				st[i,j] = st[n-1,j+int((n-1-i)/2)]
	
	s_ave[0,0][0] = s_ave_t

	for j in range(n+1):      #last term
		a_max = (b_s + s0*u*(1-(u**(n-j)))/(1-u) + s0*(u**(n-j))*d*(1-(d**j))/(1-d))/(b_n+n+1)
		a_min = (b_s + s0*d*(1-(d**j))/(1-d) + s0*(d**j)*u*(1-(u**(n-j)))/(1-u))/(b_n+n+1)
		for k in range(m+1):
			s_ave[n,j][k] = ((m-k)/m)*a_max + (k/m)*a_min
			call_value[n,j][k] = max(0, s_ave[n,j][k] - K)
		
	for i in range(n-1,-1,-1):   #backward induction
		s_ave[i,0] = (b_s + s0*u*(1-(u**i))/(1-u) )/(b_n+i+1)       #top 
		call_value[i,0] = (p*call_value[i+1,0][0] +(1-p)*call_value[i+1,1][0]) * np.exp(-r*dt)
		
		s_ave[i,i] = (b_s + s0*d*(1-(d**i))/(1-d) )/(b_n+i+1)       #bottom
		call_value[i,i] = (p*call_value[i+1,i][-1] +(1-p)*call_value[i+1,i+1][0]) * np.exp(-r*dt)
			
		for j in range(1, i):
			a_max = (b_s + s0*u*(1-(u**(i-j)))/(1-u) + s0*(u**(i-j))*d*(1-(d**j))/(1-d))/(b_n+i+1)
			a_min = (b_s + s0*d*(1-(d**j))/(1-d) + s0*(d**j)*u*(1-(u**(i-j)))/(1-u))/(b_n+i+1)
			for k in range(m+1):
				
				s_ave[i,j][k] = ((m-k)/m)*a_max + (k/m)*a_min
				
				Au = ((b_n+i+1)*s_ave[i,j][k] + st[i+1,j]) / (b_n+i+2)
				Ad = ((b_n+i+1)*s_ave[i,j][k] + st[i+1,j+1]) / (b_n+i+2)
				
				if k == 0:            #same path
					Ad = s_ave[i+1,j+1][0]
				if k == m:
					Au = s_ave[i+1,j][-1]

					
				for u_price in s_ave[i+1,j]:
					if u_price < Au:
						u_index = np.where(s_ave[i+1,j]==u_price)[0][0]
						wu = (s_ave[i+1,j][u_index-1] - Au) / (s_ave[i+1,j][u_index-1] - s_ave[i+1,j][u_index])
						cu = wu * call_value[i+1,j][u_index] + (1-wu)* call_value[i+1,j][u_index-1]
						break
					elif u_price == Au:
						u_index = np.where(s_ave[i+1,j]==u_price)[0][0]
						cu = call_value[i+1,j][u_index]
						break
								
				for d_price in s_ave[i+1,j+1]:
					if d_price < Ad:
						d_index = np.where(s_ave[i+1,j+1]==d_price)[0][0]
						wd = (s_ave[i+1,j+1][d_index-1] - Ad) / (s_ave[i+1,j+1][d_index-1] - s_ave[i+1,j+1][d_index])
						cd = wd * call_value[i+1,j+1][d_index] + (1-wd)* call_value[i+1,j+1][d_index-1]
						break
					elif d_price == Ad:
						d_index = np.where(s_ave[i+1,j+1]==d_price)[0][0]
						cd = call_value[i+1,j+1][d_index]
						break
												
				
				call_value[i,j][k] = (p*cu + (1-p)*cd) * np.exp(-r*dt)
		
	print(call_value[0,0][0])

print('european')
Euro_call(s0, r, q, sigma, T_t, t, n, s_ave_t, m, K)


def Ame_call(s0, r, q, sigma, T_t, t, n, save_t, m, K):
	dt = (T_t)/n
	u = np.exp(sigma * np.sqrt(dt))
	d = 1/u
	p = (np.exp((r-q)*(dt))-d) / (u-d)	
	b_n = (n/T_t) * t      #n before t
	b_s = s_ave_t*(b_n+1)  
	
	st = np.zeros((n+1,n+1))
	s_ave = np.zeros((n+1,n+1,m+1))
	call_value = np.zeros((n+1,n+1,m+1))
	
	st[0,0] = s0
	for j in range(0,n+1):           #n-period
		st[n,j] = s0*(u**(n-j))*(d**j)
		if n-j == j:
			st[n,j] = s0
	
	for j in range(0,n):
		st[n-1,j] = s0*(u**(n-1-j))*(d**j)
		if n-1-j == j:
			ST[n-1,j] = s0

	for i in range(1,n-1):
		for j in range(0,i+1): 
			if (n-i) % 2 == 0:
				st[i,j] = st[n,j+int((n-i)/2)]
			else:
				st[i,j] = st[n-1,j+int((n-1-i)/2)]
	
	s_ave[0,0][0] = s_ave_t

	for j in range(n+1):      #last term
		a_max = (b_s + s0*u*(1-(u**(n-j)))/(1-u) + s0*(u**(n-j))*d*(1-(d**j))/(1-d))/(b_n+n+1)
		a_min = (b_s + s0*d*(1-(d**j))/(1-d) + s0*(d**j)*u*(1-(u**(n-j)))/(1-u))/(b_n+n+1)
		for k in range(m+1):
			s_ave[n,j][k] = ((m-k)/m)*a_max + (k/m)*a_min
			call_value[n,j][k] = max(0, s_ave[n,j][k] - K)
		
	for i in range(n-1,-1,-1):   #backward induction
		s_ave[i,0] = (b_s + s0*u*(1-(u**i))/(1-u) )/(b_n+i+1)       #top 
		call_value[i,0] = max((p*call_value[i+1,0][0] +(1-p)*call_value[i+1,1][0]) * np.exp(-r*dt) , s_ave[i,0][0] -K)
		
		s_ave[i,i] = (b_s + s0*d*(1-(d**i))/(1-d) )/(b_n+i+1)       #bottom
		call_value[i,i] = max((p*call_value[i+1,i][-1] +(1-p)*call_value[i+1,i+1][0]) * np.exp(-r*dt), s_ave[i,i][0] - K)
			
		for j in range(1, i):
			a_max = (b_s + s0*u*(1-(u**(i-j)))/(1-u) + s0*(u**(i-j))*d*(1-(d**j))/(1-d))/(b_n+i+1)
			a_min = (b_s + s0*d*(1-(d**j))/(1-d) + s0*(d**j)*u*(1-(u**(i-j)))/(1-u))/(b_n+i+1)
			for k in range(m+1):
				
				s_ave[i,j][k] = ((m-k)/m)*a_max + (k/m)*a_min
				
				Au = ((b_n+i+1)*s_ave[i,j][k] + st[i+1,j]) / (b_n+i+2)
				Ad = ((b_n+i+1)*s_ave[i,j][k] + st[i+1,j+1]) / (b_n+i+2)
				
				if k == 0:            #same path
					Ad = s_ave[i+1,j+1][0]
				if k == m:
					Au = s_ave[i+1,j][-1]

					
				for u_price in s_ave[i+1,j]:
					if u_price < Au:
						u_index = np.where(s_ave[i+1,j]==u_price)[0][0]
						wu = (s_ave[i+1,j][u_index-1] - Au) / (s_ave[i+1,j][u_index-1] - s_ave[i+1,j][u_index])
						cu = wu * call_value[i+1,j][u_index] + (1-wu)* call_value[i+1,j][u_index-1]
						break
					elif u_price == Au:
						u_index = np.where(s_ave[i+1,j]==u_price)[0][0]
						cu = call_value[i+1,j][u_index]
						break
								
				for d_price in s_ave[i+1,j+1]:
					if d_price < Ad:
						d_index = np.where(s_ave[i+1,j+1]==d_price)[0][0]
						wd = (s_ave[i+1,j+1][d_index-1] - Ad) / (s_ave[i+1,j+1][d_index-1] - s_ave[i+1,j+1][d_index])
						cd = wd * call_value[i+1,j+1][d_index] + (1-wd)* call_value[i+1,j+1][d_index-1]
						break
					elif d_price == Ad:
						d_index = np.where(s_ave[i+1,j+1]==d_price)[0][0]
						cd = call_value[i+1,j+1][d_index]
						break
												
				
				call_value[i,j][k] = max((p*cu + (1-p)*cd) * np.exp(-r*dt), s_ave[i,j][k] - K)
		
	print(call_value[0,0,0])

print('american')
Ame_call(s0, r, q, sigma, T_t, t, n, s_ave_t, m, K)


def monte_carlo_ave(s0, r, q, sigma, T, t, n, s_ave_t):
	dt = T_t/n
	st = np.zeros((10000,n+1))
	st[:,0] = np.log(s0)
	call = np.zeros(20)
	b_n = (n/T_t) * t      #n before t
	b_s = s_ave_t*(b_n+1) 


	for rep in range(20):
		for i in range(n):
			e = np.random.normal(0,1,10000)
			st[:,i+1] = st[:,i] + (r-q-0.5*(sigma**2))*(dt) + sigma*((dt)**0.5)*e
		payoff = (np.exp(st[:,1:]).sum(axis=1) + b_s)/(b_n+n+1) - K
		payoff[payoff < 0] = 0
		call[rep] = np.mean(payoff)* np.exp(-r*(T_t))
	
	print(np.mean(call))
	print((np.mean(call)-np.std(call)*2), (np.mean(call)+np.std(call)*2))

print("monte carlo")
monte_carlo_ave(s0, r, q, sigma, T_t, t, n, s_ave_t)


## bonus 1
def Euro_call(s0, r, q, sigma, T_t, t, n, s_ave_t, m, K):
	dt = (T_t)/n
	u = np.exp(sigma * np.sqrt(dt))
	d = 1/u
	p = (np.exp((r-q)*(dt))-d) / (u-d)	
	b_n = (n/T_t) * t    
	b_s = s_ave_t*(b_n+1)  
	
	st = np.zeros((n+1,n+1))
	s_ave = np.zeros((n+1,n+1,m+1))
	call_value = np.zeros((n+1,n+1,m+1))
	
	st[0,0] = s0
	for j in range(0,n+1):           
		st[n,j] = s0*(u**(n-j))*(d**j)
		if n-j == j:
			st[n,j] = s0
	
	for j in range(0,n):
		st[n-1,j] = s0*(u**(n-1-j))*(d**j)
		if n-1-j == j:
			ST[n-1,j] = s0

	for i in range(1,n-1):
		for j in range(0,i+1): 
			if (n-i) % 2 == 0:
				st[i,j] = st[n,j+int((n-i)/2)]
			else:
				st[i,j] = st[n-1,j+int((n-1-i)/2)]
	
	s_ave[0,0][0] = s_ave_t

	for j in range(n+1):      #last term
		a_max = (b_s + s0*u*(1-(u**(n-j)))/(1-u) + s0*(u**(n-j))*d*(1-(d**j))/(1-d))/(b_n+n+1)
		a_min = (b_s + s0*d*(1-(d**j))/(1-d) + s0*(d**j)*u*(1-(u**(n-j)))/(1-u))/(b_n+n+1)
		for k in range(m+1):
			s_ave[n,j][k] = np.exp(((m-k)/m)*np.log(a_max)+ (k/m)*np.log(a_min))
			call_value[n,j][k] = max(0, s_ave[n,j][k] - K)
		
	for i in range(n-1,-1,-1):   #backward induction
		s_ave[i,0] = (b_s + s0*u*(1-(u**i))/(1-u) )/(b_n+i+1)     
		call_value[i,0] = (p*call_value[i+1,0][0] +(1-p)*call_value[i+1,1][0]) * np.exp(-r*dt)
		
		s_ave[i,i] = (b_s + s0*d*(1-(d**i))/(1-d) )/(b_n+i+1)     
		call_value[i,i] = (p*call_value[i+1,i][-1] +(1-p)*call_value[i+1,i+1][0]) * np.exp(-r*dt)
			
		for j in range(1, i):
			a_max = (b_s + s0*u*(1-(u**(i-j)))/(1-u) + s0*(u**(i-j))*d*(1-(d**j))/(1-d))/(b_n+i+1)
			a_min = (b_s + s0*d*(1-(d**j))/(1-d) + s0*(d**j)*u*(1-(u**(i-j)))/(1-u))/(b_n+i+1)
			for k in range(m+1):
				
				s_ave[i,j][k] = np.exp(((m-k)/m)*np.log(a_max)+ (k/m)*np.log(a_min))
				
				Au = ((b_n+i+1)*s_ave[i,j][k] + st[i+1,j]) / (b_n+i+2)
				Ad = ((b_n+i+1)*s_ave[i,j][k] + st[i+1,j+1]) / (b_n+i+2)
				
				if k == 0:            #same path
					Ad = s_ave[i+1,j+1][0]
				if k == m:
					Au = s_ave[i+1,j][-1]

					
				for u_price in s_ave[i+1,j]:
					if u_price < Au:
						u_index = np.where(s_ave[i+1,j]==u_price)[0][0]
						wu = (s_ave[i+1,j][u_index-1] - Au) / (s_ave[i+1,j][u_index-1] - s_ave[i+1,j][u_index])
						cu = wu * call_value[i+1,j][u_index] + (1-wu)* call_value[i+1,j][u_index-1]
						break
					elif u_price == Au:
						u_index = np.where(s_ave[i+1,j]==u_price)[0][0]
						cu = call_value[i+1,j][u_index]
						break
								
				for d_price in s_ave[i+1,j+1]:
					if d_price < Ad:
						d_index = np.where(s_ave[i+1,j+1]==d_price)[0][0]
						wd = (s_ave[i+1,j+1][d_index-1] - Ad) / (s_ave[i+1,j+1][d_index-1] - s_ave[i+1,j+1][d_index])
						cd = wd * call_value[i+1,j+1][d_index] + (1-wd)* call_value[i+1,j+1][d_index-1]
						break
					elif d_price == Ad:
						d_index = np.where(s_ave[i+1,j+1]==d_price)[0][0]
						cd = call_value[i+1,j+1][d_index]
						break
												
				
				call_value[i,j][k] = (p*cu + (1-p)*cd) * np.exp(-r*dt)
		
	print(call_value[0,0][0])

##print('bonus_european')
#Euro_call(s0, r, q, sigma, T_t, t, n, s_ave_t, m, K)


def Ame_call(s0, r, q, sigma, T_t, t, n, save_t, m, K):
	dt = (T_t)/n
	u = np.exp(sigma * np.sqrt(dt))
	d = 1/u
	p = (np.exp((r-q)*(dt))-d) / (u-d)	
	b_n = (n/T_t) * t      #n before t
	b_s = s_ave_t*(b_n+1)  
	
	st = np.zeros((n+1,n+1))
	s_ave = np.zeros((n+1,n+1,m+1))
	call_value = np.zeros((n+1,n+1,m+1))
	
	st[0,0] = s0
	for j in range(0,n+1):           #n-period
		st[n,j] = s0*(u**(n-j))*(d**j)
		if n-j == j:
			st[n,j] = s0
	
	for j in range(0,n):
		st[n-1,j] = s0*(u**(n-1-j))*(d**j)
		if n-1-j == j:
			ST[n-1,j] = s0

	for i in range(1,n-1):
		for j in range(0,i+1): 
			if (n-i) % 2 == 0:
				st[i,j] = st[n,j+int((n-i)/2)]
			else:
				st[i,j] = st[n-1,j+int((n-1-i)/2)]
	
	s_ave[0,0][0] = s_ave_t

	for j in range(n+1):      #last term
		a_max = (b_s + s0*u*(1-(u**(n-j)))/(1-u) + s0*(u**(n-j))*d*(1-(d**j))/(1-d))/(b_n+n+1)
		a_min = (b_s + s0*d*(1-(d**j))/(1-d) + s0*(d**j)*u*(1-(u**(n-j)))/(1-u))/(b_n+n+1)
		for k in range(m+1):
			s_ave[n,j][k] = np.exp(((m-k)/m)*np.log(a_max)+ (k/m)*np.log(a_min))
			call_value[n,j][k] = max(0, s_ave[n,j][k] - K)
		
	for i in range(n-1,-1,-1):   #backward induction
		s_ave[i,0] = (b_s + s0*u*(1-(u**i))/(1-u) )/(b_n+i+1)       #top 
		call_value[i,0] = max((p*call_value[i+1,0][0] +(1-p)*call_value[i+1,1][0]) * np.exp(-r*dt) , s_ave[i,0][0] -K)
		
		s_ave[i,i] = (b_s + s0*d*(1-(d**i))/(1-d) )/(b_n+i+1)       #bottom
		call_value[i,i] = max((p*call_value[i+1,i][-1] +(1-p)*call_value[i+1,i+1][0]) * np.exp(-r*dt), s_ave[i,i][0] - K)
			
		for j in range(1, i):
			a_max = (b_s + s0*u*(1-(u**(i-j)))/(1-u) + s0*(u**(i-j))*d*(1-(d**j))/(1-d))/(b_n+i+1)
			a_min = (b_s + s0*d*(1-(d**j))/(1-d) + s0*(d**j)*u*(1-(u**(i-j)))/(1-u))/(b_n+i+1)
			for k in range(m+1):
				
				s_ave[i,j][k] = np.exp(((m-k)/m)*np.log(a_max)+ (k/m)*np.log(a_min))
				
				Au = ((b_n+i+1)*s_ave[i,j][k] + st[i+1,j]) / (b_n+i+2)
				Ad = ((b_n+i+1)*s_ave[i,j][k] + st[i+1,j+1]) / (b_n+i+2)
				
				if k == 0:            #same path
					Ad = s_ave[i+1,j+1][0]
				if k == m:
					Au = s_ave[i+1,j][-1]

					
				for u_price in s_ave[i+1,j]:
					if u_price < Au:
						u_index = np.where(s_ave[i+1,j]==u_price)[0][0]
						wu = (s_ave[i+1,j][u_index-1] - Au) / (s_ave[i+1,j][u_index-1] - s_ave[i+1,j][u_index])
						cu = wu * call_value[i+1,j][u_index] + (1-wu)* call_value[i+1,j][u_index-1]
						break
					elif u_price == Au:
						u_index = np.where(s_ave[i+1,j]==u_price)[0][0]
						cu = call_value[i+1,j][u_index]
						break
								
				for d_price in s_ave[i+1,j+1]:
					if d_price < Ad:
						d_index = np.where(s_ave[i+1,j+1]==d_price)[0][0]
						wd = (s_ave[i+1,j+1][d_index-1] - Ad) / (s_ave[i+1,j+1][d_index-1] - s_ave[i+1,j+1][d_index])
						cd = wd * call_value[i+1,j+1][d_index] + (1-wd)* call_value[i+1,j+1][d_index-1]
						break
					elif d_price == Ad:
						d_index = np.where(s_ave[i+1,j+1]==d_price)[0][0]
						cd = call_value[i+1,j+1][d_index]
						break
												
				
				call_value[i,j][k] = max((p*cu + (1-p)*cd) * np.exp(-r*dt), s_ave[i,j][k] - K)
		
	print(call_value[0,0][0])

#  print('bonus_american')
#Ame_call(s0, r, q, sigma, T_t, t, n, s_ave_t, m, K)








