d=3
X=[0,1,0,0.3,0,1,0,0,1,0]
Y=[0,0,1,0.4,0,0,1,0,0,1]
for i in range(0,10):
	x=X[i]
	y=Y[i]
	for d in range(3,-1,-1):
		for px in range(d,-1,-1):
			if i <= 4:
				print(pow(x,px)*pow(y,d-px), end =" ")
			elif i <= 7:
				if px==0:print(0, end =" ")
				else:print(px*pow(x,px-1)*pow(y,d-px), end =" ")
			else:
				if d-px==0:print(0, end =" ")
				else:print((d-px)*pow(x,px)*pow(y,d-px-1), end =" ")
			print("&",end=" ")
	print('\\\\')

