#This code implements Algorithm 4.2 of ``An enriched count of the bitangents to a smooth plane quartic curve"
#by Hannah Larson and Isabel Vogt

#Given a real polynomial f in x, y, z defining a smooth plane quartic, all_counts(bitangents_list(f)) returns 
#a list of all possible signed counts of bitangents to V(f) obtained by varying the line at infinity. 

T.<x,y,z> = PolynomialRing(QQ)

#All computations are done using floating point;
#It may be necessary to increase the precision when the coefficents are large.

precis = 10000
threshold = 100
RRR = RealField(precis)

#The Trott curve:
f = 12^2 *(x^4 + y^4) - 15^2*(x^2 + y^2)*z^2 + 350*x^2*y^2 + 81*z^4

print "f is", f

#INPUT: floating point real polynomial P
#OUTPUT: list of floating point real roots
def FindRealRoots(P):

	Roots =[]
	AllRoots = P.complex_roots()
	for w in AllRoots:
		if abs(w.imag()) < 2^(-threshold):
			found = False
			for u in Roots:
				if abs(u - w) < 2^(-threshold):
					found = True
			if not found:
				Roots.append(w.real())
	return Roots


#INPUT: a homogeneous quartic in QQ[x,y,z]
#OUTPUT: a list [real bitangent, [affine points of bitangency], type] if the points are real, or [real bitangent, FALSE, 1] if they are complex conjugate
#ASSUMING FOR NOW THAT NONE OF THE BITANGENTS OCCUR ALONG THE LINE z=0

def bitangents_list(f):
	
	T.<x,y,z> = PolynomialRing(QQ)

	dfx = f.derivative(x).subs(z = 1)
	dfy = f.derivative(y).subs(z = 1)
	dfz = f.derivative(z).subs(z = 1)

	#FIRST FIND THE LINES OF THE FORM Y = R*X + S*Z
	S.<r,s,u> = PolynomialRing(QQ)
	R.<x,y,z> = PolynomialRing(S)
	f0 = f.base_extend(R)
	
	f0 = f0.subs(y = r*x + s*z)
	f0 = f0.subs(z = 1)
	a = f0.coefficient(x^4).constant_coefficient()
	b = f0.coefficient(x^3).constant_coefficient()
	c = f0.coefficient(x^2).constant_coefficient()
	d = f0.coefficient(x).constant_coefficient()
	e = f0.constant_coefficient()

	if a == 0:
		raise ValueError("bitangent points along infinity")
	alpha = b/(2*a)
	beta = c/(2*a) - b^2/(8*a^2)
	#The remaining two equations for d and e are:
	
	#2*a*alpha*beta = d, which gives
	eq1 = (8*a^2*d - 4*a*b*c + b^3)

	#a*beta^2 = e, which gives
	eq2 = (64*a^3*e - 16*a^2*c^2 + 8*a*b^2*c - b^4)	

	aResb = a.resultant(b,r)
	eq1Reseq2 = eq1.resultant(eq2, r)

	
	SS.<s> = PolynomialRing(QQ)

	aResb = SS(aResb)
	eq1Reseq2 = SS(eq1Reseq2)

	Res = SS(eq1Reseq2/aResb^8)
	Res = Res/Res[Res.degree()]

	
	
	Res0 = Res.base_extend(RRR)	
	RRoots_s = FindRealRoots(Res0)

	RBit = []
	for s0 in RRoots_s:
		eq1r = eq1.subs(s = s0)
		eq1r = eq1r.univariate_polynomial()
		if eq1r == 0:
			eq2r = eq2.subs(s=s0)
			eq2r = eq2r.univariate_polynomial()
			
			for r0 in FindRealRoots(eq2r):
				RBit.append([s0, r0])
		else:
			RRoots_r = FindRealRoots(eq1r)
			for r0 in RRoots_r:
				if abs(eq2.subs(s = s0, r = r0)) < 2^(-100):
					RBit.append([s0, r0])
	
	
	#Now we want to find the bitangent line and bitangency points
	BPoints = []
	for ll in RBit:
		ss = ll[0]
		rr = ll[1]
		aalpha = alpha.subs(r = rr, s = ss)
		bbeta = beta.subs(r = rr, s = ss)
		if aalpha^2 - 4*bbeta < 2^(-threshold):
			T.<x,y,z>=PolynomialRing(RRR)
			L = -rr*x + y - ss*z
			L = T(sum([T.gens()[i] * L.coefficient(L.parent().gens()[i]).constant_coefficient() for i in xrange(3)]))
			BPoints.append([L, False, 1])
		
		elif aalpha^2 - 4*bbeta >= 2^(-threshold):
			pts_equation = (x^2 + aalpha*x + bbeta).univariate_polynomial()
			pts_list = []
			for x0 in FindRealRoots(pts_equation):
				pts_list.append((x0, rr*x0 + ss))
	                T.<x,y,z>=PolynomialRing(RRR)
			L = -rr*x + y - ss*z
	                L = T(sum([T.gens()[i] * L.coefficient(L.parent().gens()[i]).constant_coefficient() for i in xrange(3)]))
			dfL = L.coefficient(x)*dfx + L.coefficient(y)*dfy + L.coefficient(z)*dfz

        		if dfL.subs(x = pts_list[0][0], y = pts_list[0][1], z=1)*dfL.subs(x = pts_list[1][0], y = pts_list[1][1], z = 1) > 0:
                		t = 1
        		else:
                		t = -1
			BPoints.append([L, pts_list, t])
	
	#NOW LET'S FIND THE LINES WHERE THE COEFF OF Y IS ZERO
	S.<r,s,u> = PolynomialRing(QQ)
	R.<x,y,z> = PolynomialRing(S)
	f1 = f.base_extend(R)
	
	f1 = f1.subs(x = u*z)
	f1 = f1.subs(z = 1)
	a1 = f1.coefficient(y^4).constant_coefficient()
	b1 = f1.coefficient(y^3).constant_coefficient()
	c1 = f1.coefficient(y^2).constant_coefficient()
	d1 = f1.coefficient(y).constant_coefficient()
	e1 = f1.constant_coefficient()

	if a1 == 0:
		raise ValueError("bitangent along infinity")
	
	#In order for a1*y^4 + b1*y^3 + c1*y^2 + d1*y + e1 = a1(y^2 + alpha1*y + beta)^2, we must have:
	alpha1 = b1/(2*a1)
	beta1 = c1/(2*a1) - b1^2/(8*a1^2)
	#The reaming two equations for d and e are:
	
	#2*a1*alpha1*beta1 = d1, which gives
	eq11 = (8*a1^2*d1 - 4*a1*b1*c1 + b1^3)
	eq11 = eq11.univariate_polynomial()
	
	#a1*beta1^2 = e1, which gives
	eq21 = (64*a1^3*e1 - 16*a1^2*c1^2 + 8*a1*b1^2*c1 - b1^4)
	eq21 = eq21.univariate_polynomial()
	
	RBit2 = []
	
	eq110 = eq11.base_extend(RRR)
	if eq110 == 0:
		eq210 = eq21.base_extend(RRR)
		for uu in FindRealRoots(eq210):
			RBit2.append(uu)
	else:
		RRoots_u = FindRealRoots(eq110)
		for uu in RRoots_u:
			if abs(eq21.subs(u = uu)) < 2^(-100):
				RBit2.append(uu)
	
	for uu in RBit2:
	        aalpha1 = alpha1.subs(u = uu)
	        bbeta1 = beta1.subs(u = uu)
	        if aalpha1^2 - 4*bbeta1 < 2^(-threshold):
	                T.<x,y,z>=PolynomialRing(RRR)
	                L = x - uu*z
	                L = T(sum([T.gens()[i] * L.coefficient(L.parent().gens()[i]).constant_coefficient() for i in xrange(3)]))
	                BPoints.append([L, False, 1])
	
	        elif aalpha1^2 - 4*bbeta1 >= 2^(-threshold):
	                pts_equation = (y^2 + aalpha1*y + bbeta1).univariate_polynomial()
	                pts_list = []
	                for y0 in FindRealRoots(pts_equation):
	                        pts_list.append((uu, y0))
	                T.<x,y,z>=PolynomialRing(RRR)
	                L = x - uu*z
	                L = T(sum([T.gens()[i] * L.coefficient(L.parent().gens()[i]).constant_coefficient() for i in xrange(3)]))
	                dfL = L.coefficient(x)*dfx + L.coefficient(y)*dfy + L.coefficient(z)*dfz
                        
                        if dfL.subs(x = pts_list[0][0], y = pts_list[0][1], z=1)*dfL.subs(x = pts_list[1][0], y = pts_list[1][1], z = 1) > 0:
                                t = 1
                        else:
                                t = -1
			BPoints.append([L, pts_list,t])

	print "The number of real bitangents is", len(BPoints)

	return BPoints


#INPUT: the output of bitangents_list, e.g., a list of tuples [real bitangent line, affine points of bitangency, type relative to z=0]
#OUTPUT: signed count relative to L_infty = V(z)
def signed_count(realBlist):	

        cc_lines = 28 - len(realBlist)
        pcounter = cc_lines/2
        ncounter = cc_lines/2

        for L in realBlist:
                if L[2] == 1:
                        pcounter += 1
                elif L[2] == -1:
                        ncounter += 1

        print "The signed count is ", pcounter, "<1> + ", ncounter, "<-1>"
        return pcounter - ncounter

#INPUT: the output of bitangents_list, e.g., a list of tuples [real bitangent line, affine points of bitangency, type relative to z=0]
#OUTPUT: a set of all possible signed counts for varying lines at infinity, following Algorithm 4.2

def all_counts(realBlist):
	
	#The signed count when L_infty is the line z=0
	original = signed_count(realBlist)


	#The two points of bitangency partition the bitangent line into two real components 
	#(one of which may be empty if it is a hyperflex)
	#The segment that does not intersect the line z=0 (the original L_infty) is called a GRATE
	#The signed count relative to an aribitary L_infty can only differ from the signed count
	#relative to z=0 if L_infty intersects one of the grates.

	#We first make the list of grates as [type rel z=0, (x1,y1), (x2, y2)]
	#to indicate that the grate is the line segment between (x1,y1) and (x2, y2) in
	#the affine plane with coordinates x and y

	grate = []
	for bline in realBlist:
		if not bline[1] == False:
			grate.append([bline[2], bline[1][0], bline[1][1]])

	#We now make a list of x-coordinates of grate points, together with the 
	#the effect on the type if a line at infinty flips that grate

	xlist = [] #store pairs, (x coord, effect)
	for [t, (x1,y1), (x2, y2)] in grate:
        	if abs(x1 - x2) > 2^(-threshold):
                	xlist.append((min(x1,x2), -2*t))
                	xlist.append((max(x1,x2), 2*t))

	#This makes the possible type for the pencil of lines through the point at infinity (e.g.,
	#vertical in the affine plane.

	xlist.sort()
	poss = [original]
	count = original
	for (x, e) in xlist:
        	count = count + e
        	poss.append(count)
	print set(poss)

	#We now sample the finitely many other pencils (corresponding to red lines in Figure 1)

	grate_points = []
	for [t, (x1,y1), (x2, y2)] in grate:
        	grate_points.append((x1,y1))
        	grate_points.append((x2,y2))
	
	#First store the slopes corresponding to black dashed lines as in Figure 1 in Algorithm 4.2
	
	slope_list = []
	for i in [0..len(grate_points)-1]:
        	(x1,y1) = grate_points[i]
        	for j in [i+1..len(grate_points)-1]:
                	(x2,y2) = grate_points[j]
                	slope_list.append((y2-y1)/(x2-x1))

	slope_list.sort()
	
	#The midpoints of slope_list will give slopes corresponding to the red lines in Figure 1

	midpoints = []
	for i in [0..len(slope_list)-2]:
        	midpoints.append((slope_list[i] + slope_list[i+1])/2)

	#Now check along each red test pencil

	for s in midpoints:
        	count = original
        	xlist = []
		this_poss = []
        	for [t, (x1,y1), (x2, y2)] in grate:
                	x1 = (s*x1 - y1)/(1/s + s)
                	x2 = (s*x2 - y2)/(1/s + s)
                	if abs(x1 - x2) > 2^(-threshold):
				xlist.append((min(x1,x2), -2*t))
                		xlist.append((max(x1,x2), 2*t))
        	xlist.sort()
        	for (x, e) in xlist:
                	count = count + e
                	poss.append(count)
			
		

	return set(poss)


