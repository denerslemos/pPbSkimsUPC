int get_Ntrkoff(int size, float *eta, float *pt, int *charge, bool *hp, float *pterr, float *dcaxy, float *dcaxyerr,  float *dcaz, float *dcazerr){
	int Ntrk_off = 0;
	for(int ii=0; ii<size; ii++){ 
		if(pt[ii] <= 0.4) continue;
		if(fabs(eta[ii]) >= 2.4) continue; 
		if(fabs(charge[ii]) == 0)continue;
		if(hp[ii] != 1) continue;
		if(fabs(pterr[ii]/pt[ii]) >= 0.1) continue;
		if(fabs(dcaxy[ii]/dcaxyerr[ii]) >= 3.0) continue;
		if(fabs(dcaz[ii]/dcazerr[ii]) >= 3.0) continue;
		Ntrk_off=Ntrk_off+1;
	}
	return Ntrk_off;
}
