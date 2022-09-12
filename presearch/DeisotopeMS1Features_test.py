from numpy import *
import copy

def DeisotopeMS1Features(monomz,monochg,monointen,monoppi,mz,inten,unitdiff,ptol,all_tIPV):
	# FilterByModel
	
	mz_len=len(monomz)
	Dscore=array([0.0]*mz_len)
	delta_Dscore=array([0.0]*mz_len)
	
	#Dscore
	for ino in range(mz_len):
		c_pmz=monomz[ino]
		c_chg=monochg[ino]
		
		c_mz = [c_pmz-unitdiff/c_chg,c_pmz,c_pmz+unitdiff/c_chg,c_pmz+2*unitdiff/c_chg,c_pmz+3*unitdiff/c_chg,c_pmz+4*unitdiff/c_chg,c_pmz+5*unitdiff/c_chg]
		[e_mz,e_inten] = AlignMS1(mz,inten,c_mz,ptol)
		if abs(c_pmz-669.39)<=0.1 and c_chg==3:
			print('---e_inten---')
			for i in range(0,len(e_mz)):
				print('%.4f %.0f' % (e_mz[i],e_inten[i]))
		
		
		Dscore[ino]=get_sim(c_pmz,c_chg,e_inten,all_tIPV)
		monointen[ino]=sum(e_inten)
	print('---monomz---')
	for i in range(0,len(monomz)):
		print('%.4f %d %.3f' % (monomz[i],monochg[i],Dscore[i]))
	#delta_Dscore
	max_Dscore=max(Dscore)+1e-10
	for ino in range(mz_len):
		delta_Dscore[ino]=(max_Dscore-Dscore[ino])/max_Dscore
	
	idx=nonzero(delta_Dscore<0.1)[0]
	if len(idx)>0:
		monomz=array(monomz)
		monochg=array(monochg)
		monointen=array(monointen)
		
		monomz=monomz[idx]
		monochg=monochg[idx]
		monointen=monointen[idx]
		sum_inten=sum(monointen)+0.1
		monoppi=monointen/sum_inten
		
		monomz=monomz.tolist()
		monochg=monochg.tolist()
		monointen=monointen.tolist()
		monoppi=monoppi.tolist()
	
	return [monomz,monochg,monointen,monoppi]

def AlignMS1(mz,inten,mz_n,ptol):
	
	mz_t=[]
	inten_t=[]
	
	for ino in range(0,len(mz_n)):
		c_ptol=ptol*mz_n[ino]*1e-6
		left=mz_n[ino]-c_ptol
		right=mz_n[ino]+c_ptol
		pos=nonzero((mz>=left) & (mz<=right))[0]
		if len(pos)==0:
			mz_t.append(mz_n[ino])
			inten_t.append(0.0)
		else:
			mz_t.append(sum(mz[pos]*inten[pos])/sum(inten[pos]))
			inten_t.append(max(inten[pos]))
	
	mz_t=array(mz_t)
	inten_t=array(inten_t)
	
	return [mz_t,inten_t]

def get_sim(c_pmz,c_chg,e_inten,all_tIPV):
	#
	pmass=1.007276
	M = c_pmz*c_chg-c_chg*pmass
	
	if len(all_tIPV)==1:
		TMT_data=0
	else:
		TMT_data=1
	
	if TMT_data==0:
		all_tIPV_0 = all_tIPV["0"]
		tIPV_0=all_tIPV_0[min(len(all_tIPV_0)-1,int(M)-1)]
		
		Dscore=get_Dscore(e_inten,tIPV_0)
	else:
		all_tIPV_0 = all_tIPV["0"]
		tIPV_0=all_tIPV_0[min(len(all_tIPV_0)-1,int(M)-1)]
		all_tIPV_1 = all_tIPV["1"]
		tIPV_1=all_tIPV_1[min(len(all_tIPV_1)-1,int(M)-1)]
		all_tIPV_2 = all_tIPV["2"]
		tIPV_2=all_tIPV_2[min(len(all_tIPV_2)-1,int(M)-1)]
		all_tIPV_3 = all_tIPV["3"]
		tIPV_3=all_tIPV_3[min(len(all_tIPV_3)-1,int(M)-1)]
		all_tIPV_4 = all_tIPV["4"]
		tIPV_4=all_tIPV_4[min(len(all_tIPV_4)-1,int(M)-1)]
		
		Dscore0=get_Dscore(e_inten,tIPV_0)
		Dscore1=get_Dscore(e_inten,tIPV_1)
		Dscore2=get_Dscore(e_inten,tIPV_2)
		Dscore3=get_Dscore(e_inten,tIPV_3)
		Dscore4=get_Dscore(e_inten,tIPV_4)
		
		Dscore=max([Dscore0,Dscore1,Dscore2,Dscore3,Dscore4])
		if abs(c_pmz-669.39)<=0.1 and c_chg==3:
			print('---t_inten0---')
			for i in range(0,len(tIPV_0)):
				print(' %.0f' % (tIPV_0[i]))
			print('---t_inten1---')
			for i in range(0,len(tIPV_1)):
				print(' %.0f' % (tIPV_1[i]))
			print('---t_inten2---')
			for i in range(0,len(tIPV_2)):
				print(' %.0f' % (tIPV_2[i]))
			print('---t_inten3---')
			for i in range(0,len(tIPV_3)):
				print(' %.0f' % (tIPV_3[i]))
			print('---t_inten4---')
			for i in range(0,len(tIPV_4)):
				print(' %.0f' % (tIPV_4[i]))
			print('---Dscore---')
			print('%.3f %.3f %.3f %.3f %.3f, max:%.3f' % (Dscore0,Dscore1,Dscore2,Dscore3,Dscore4,Dscore))
	
	return Dscore

def get_Dscore(e_inten,tIPV):
	
	nlim=len(tIPV)
	t_inten=array([0.0]*nlim)
	for tno in range(0,nlim):
		t_inten[tno] = tIPV[tno]
	
	# cut short clusters by relative intensity
	for tno in range(3,nlim):
		if t_inten[tno]<10:
			tm=tno
			break
	if t_inten[nlim-1]>=10:
		tm=nlim
	e_inten=e_inten[0:tm]
	t_inten=t_inten[0:tm]
	
	Dscore0_global=get_similarity(e_inten,t_inten)
	Dscore0_local=get_similarity(e_inten[0:3],t_inten[0:3])
	Dscore0=0.5*Dscore0_global+0.5*Dscore0_local
	
	return Dscore0

def get_similarity(e_intens,t_intens):
	
	if sqrt( sum(e_intens*e_intens)*sum(t_intens*t_intens) )==0:
		e_sim=0
	else:
		e_sim = sum(e_intens*t_intens)/sqrt( sum(e_intens*e_intens)*sum(t_intens*t_intens) )
	
	return e_sim
