/** @file transfer.c Documented transfer module.
 *
 * Julien Lesgourgues, 28.07.2013
 *
 * This module has two purposes:
 *
 * - at the beginning, to compute the transfer functions \f$
 *   \Delta_l^{X} (q) \f$, and store them in tables used for
 *   interpolation in other modules.
 *
 * - at any time in the code, to evaluate the transfer functions (for
 *   a given mode, initial condition, type and multipole l) at any
 *   wavenumber q (by interpolating within the interpolation table).
 *
 * Hence the following functions can be called from other modules:
 *
 * -# transfer_init() at the beginning (but after perturb_init()
 *    and bessel_init())
 *
 * -# transfer_functions_at_q() at any later time
 *
 * -# transfer_free() at the end, when no more calls to
 *    transfer_functions_at_q() are needed
 *
 * Note that in the standard implementation of CLASS, only the pre-computed
 * values of the transfer functions are used, no interpolation is necessary;
 * hence the routine transfer_functions_at_q() is actually never called.
 */

#include "transfer.h"

#define AERF1 0.0705230784
#define AERF2 0.0422820123
#define AERF3 0.0092705272
#define AERF4 0.0001520143
#define AERF5 0.0002765672
#define AERF6 0.0000430638

double erf_mine(double x)
{
  if(x==0.)
    return 0;
  else {
    double y,pol;
    y=fabs(x);
    pol=(((((AERF6*y+AERF5)*y+AERF4)*y+AERF3)*y+AERF2)*y+AERF1)*y+1;
    pol=pol*pol; //2
    pol=pol*pol; //4
    pol=pol*pol; //8
    pol=pol*pol; //16
    
    if(x<0)
      return 1./pol-1;
    else
      return 1-1./pol;
  }
}

void free_selection_functions_dam(
				  struct perturbs * ppt,
				  struct transfers * ptr
				  )
{
  int ii;
 
  if(ppt->has_weak_lensing) {
    for(ii=0;ii<ptr->selection_num_wl;ii++) {
      if(ppt->has_intrinsic_alignment) {
	free(ptr->glb_selection_0_wl[ii]);
	free(ptr->glb_selection_0_tau0_minus_tau_wl[ii]);
      }
      if(ppt->has_lensing_shear) {
	free(ptr->glb_selection_lns_wl[ii]);
	free(ptr->glb_selection_lns_tau0_minus_tau_wl[ii]);
      }
    }
    if(ppt->has_intrinsic_alignment) {
      free(ptr->glb_selection_0_wl);
      free(ptr->glb_selection_0_tau0_minus_tau_wl);
      free(ptr->glb_selection_0_size_wl);
    }
    if(ppt->has_lensing_shear) {
      free(ptr->glb_selection_lns_wl);
      free(ptr->glb_selection_lns_tau0_minus_tau_wl);
      free(ptr->glb_selection_lns_size_wl);
    }
  }

  if(ppt->has_cl_number_count) {
    for(ii=0;ii<ptr->selection_num_nc;ii++) {
      if((ppt->has_nc_density) ||
	 (ppt->has_nc_rsd1) || (ppt->has_nc_rsd2) || (ppt->has_nc_rsd3) ||
	 (ppt->has_nc_gr1) || (ppt->has_nc_gr2) || (ppt->has_nc_gr3)) {
	free(ptr->glb_selection_0[ii]);
	free(ptr->glb_selection_0_tau0_minus_tau[ii]);
      }
      
      if(ppt->has_nc_lens)
	free(ptr->glb_selection_lns_lns[ii]);
      if(ppt->has_nc_gr4)
	free(ptr->glb_selection_lns_gr4[ii]);
      if(ppt->has_nc_gr5)
	free(ptr->glb_selection_lns_gr5[ii]);
      if((ppt->has_nc_lens) || (ppt->has_nc_gr4) || (ppt->has_nc_gr5)) {
	free(ptr->glb_selection_lns_tau0_minus_tau_nc[ii]);
      }
    }

    if((ppt->has_nc_density) ||
       (ppt->has_nc_rsd1) || (ppt->has_nc_rsd2) || (ppt->has_nc_rsd3) ||
       (ppt->has_nc_gr1) || (ppt->has_nc_gr2) || (ppt->has_nc_gr3)) {
      free(ptr->glb_selection_0);
      free(ptr->glb_selection_0_tau0_minus_tau);
      free(ptr->glb_selection_0_size);
    }
    
    if(ppt->has_nc_lens)
      free(ptr->glb_selection_lns_lns);
    if(ppt->has_nc_gr4)
      free(ptr->glb_selection_lns_gr4);
    if(ppt->has_nc_gr5)
      free(ptr->glb_selection_lns_gr5);
    if((ppt->has_nc_lens) || (ppt->has_nc_gr4) || (ppt->has_nc_gr5)) {
      free(ptr->glb_selection_lns_size_nc);
      free(ptr->glb_selection_lns_tau0_minus_tau_nc);
    }
  }
}

int init_selection_functions_lensing_wl_dam(
				 struct precision * ppr,
				 struct background * pba,
				 struct perturbs * ppt,
				 struct transfers * ptr
				 ) {
  double *pvecback,*selection,*tau0_minus_tau_lns_sources,*w_trapz_lns_sources;
  int ibin,tau_sources_size;
  double tau0=pba->conformal_age;

  class_alloc(ptr->glb_selection_lns_wl,ptr->selection_num_wl*sizeof(double *),ptr->error_message);
  class_alloc(ptr->glb_selection_lns_tau0_minus_tau_wl,ptr->selection_num_wl*sizeof(double *),ptr->error_message);
  class_alloc(ptr->glb_selection_lns_size_wl,ptr->selection_num_wl*sizeof(int),ptr->error_message);
  
  tau_sources_size=ppr->selection_sampling;
  class_alloc(pvecback,pba->bg_size*sizeof(double),ptr->error_message);
  class_alloc(selection,tau_sources_size*sizeof(double),ptr->error_message);
  class_alloc(tau0_minus_tau_lns_sources,tau_sources_size*sizeof(double),
	      ptr->error_message);
  class_alloc(w_trapz_lns_sources,tau_sources_size*sizeof(double),
	      ptr->error_message);

  for(ibin=0;ibin<ptr->selection_num_wl;ibin++) {
    double tau_min,tau_mean,tau_max;
    int tau_size,l_limber,index_tau;

    //Compute tau limits and size
    class_call(transfer_selection_times_wl(ppr,pba,ppt,ptr,ibin,
					   &tau_min,&tau_mean,&tau_max),
	       ptr->error_message,ptr->error_message);

    tau_size=ppr->selection_sampling;
    l_limber=ppr->l_switch_limber_for_cl_lensing;
    tau_size=MAX(tau_size,(int)((tau0-tau_min)/((tau0-tau_mean)/2./l_limber))*
		 ppr->selection_sampling_bessel);

    //Allocate and assign sizes
    ptr->glb_selection_lns_size_wl[ibin]=tau_size;
    class_alloc(ptr->glb_selection_lns_wl[ibin],tau_size*sizeof(double),ptr->error_message);
    class_alloc(ptr->glb_selection_lns_tau0_minus_tau_wl[ibin],tau_size*sizeof(double),ptr->error_message);

    //Get sampling values of tau_lensing
    class_call(transfer_selection_sampling_wl(ppr,pba,ppt,ptr,ibin,
					      tau0_minus_tau_lns_sources,tau_sources_size),
	       ptr->error_message,ptr->error_message);
    class_call(array_trapezoidal_mweights(tau0_minus_tau_lns_sources,tau_sources_size,
					  w_trapz_lns_sources,ptr->error_message),
	       ptr->error_message,ptr->error_message);
    class_call(transfer_selection_compute_wl(ppr,pba,ppt,ptr,selection,
					     tau0_minus_tau_lns_sources,
					     w_trapz_lns_sources,tau_sources_size,
					     pvecback,tau0,ibin),
	       ptr->error_message,ptr->error_message);

    //Get sampling values of tau
    class_call(transfer_lensing_sampling_wl(ppr,pba,ppt,ptr,ibin,tau0,
					    ptr->glb_selection_lns_tau0_minus_tau_wl[ibin],
					    ptr->glb_selection_lns_size_wl[ibin]),
	       ptr->error_message,ptr->error_message);

    for(index_tau=0;index_tau<tau_size;index_tau++) {
      double rescaling_wl=0;

      if(index_tau==tau_size-1)
	rescaling_wl=0;
      else {
	int index_tau_sources;

	rescaling_wl=0;
	
	for(index_tau_sources=0;index_tau_sources<tau_sources_size;index_tau_sources++) {
	  if((tau0_minus_tau_lns_sources[index_tau_sources]>0.) &&
	     (tau0_minus_tau_lns_sources[index_tau_sources]-
	      ptr->glb_selection_lns_tau0_minus_tau_wl[ibin][index_tau]>0.)) {
	    rescaling_wl+=0.5*(ptr->glb_selection_lns_tau0_minus_tau_wl[ibin][index_tau]-
			       tau0_minus_tau_lns_sources[index_tau_sources])
	      /ptr->glb_selection_lns_tau0_minus_tau_wl[ibin][index_tau]
	      /tau0_minus_tau_lns_sources[index_tau_sources]
	      *selection[index_tau_sources]
	      *w_trapz_lns_sources[index_tau_sources];
	  }
	}
      }

      ptr->glb_selection_lns_wl[ibin][index_tau]=rescaling_wl;
    }
  }

  free(pvecback);
  free(selection);
  free(tau0_minus_tau_lns_sources);
  free(w_trapz_lns_sources);

  return _SUCCESS_;
}

int init_selection_functions_lensing_nc_dam(
				 struct precision * ppr,
				 struct background * pba,
				 struct perturbs * ppt,
				 struct transfers * ptr
				 ) {
  double *pvecback,*selection,*tau0_minus_tau_lns_sources,*w_trapz_lns_sources;
  int ibin,tau_sources_size,last_index=0;
  double tau0=pba->conformal_age;

  if(ppt->has_nc_lens)
    class_alloc(ptr->glb_selection_lns_lns,ptr->selection_num_nc*sizeof(double *),ptr->error_message);
  if(ppt->has_nc_gr4)
    class_alloc(ptr->glb_selection_lns_gr4,ptr->selection_num_nc*sizeof(double *),ptr->error_message);
  if(ppt->has_nc_gr5)
    class_alloc(ptr->glb_selection_lns_gr5,ptr->selection_num_nc*sizeof(double *),ptr->error_message);
  class_alloc(ptr->glb_selection_lns_tau0_minus_tau_nc,ptr->selection_num_nc*sizeof(double *),ptr->error_message);
  class_alloc(ptr->glb_selection_lns_size_nc,ptr->selection_num_nc*sizeof(int),ptr->error_message);
  
  tau_sources_size=ppr->selection_sampling;
  class_alloc(pvecback,pba->bg_size*sizeof(double),ptr->error_message);
  class_alloc(selection,tau_sources_size*sizeof(double),ptr->error_message);
  class_alloc(tau0_minus_tau_lns_sources,tau_sources_size*sizeof(double),
	      ptr->error_message);
  class_alloc(w_trapz_lns_sources,tau_sources_size*sizeof(double),
	      ptr->error_message);

  for(ibin=0;ibin<ptr->selection_num_nc;ibin++) {
    double tau_min,tau_mean,tau_max;
    int tau_size,l_limber,index_tau;
    double *array_bg=NULL;
    short itr=ptr->selection_tracer_nc[ibin];

    //Compute tau limits and size
    class_call(transfer_selection_times_nc(ppr,pba,ppt,ptr,ibin,
					&tau_min,&tau_mean,&tau_max),
	       ptr->error_message,ptr->error_message);

    tau_size=ppr->selection_sampling;
    l_limber=ppr->l_switch_limber_for_cl_density;
    tau_size=MAX(tau_size,(int)((tau0-tau_min)/((tau0-tau_mean)/2./l_limber))*
		 ppr->selection_sampling_bessel);

    //Allocate and assign sizes
    ptr->glb_selection_lns_size_nc[ibin]=tau_size;
    if(ppt->has_nc_lens)
      class_alloc(ptr->glb_selection_lns_lns[ibin],tau_size*sizeof(double),ptr->error_message);
    if(ppt->has_nc_gr4)
      class_alloc(ptr->glb_selection_lns_gr4[ibin],tau_size*sizeof(double),ptr->error_message);
    if(ppt->has_nc_gr5)
      class_alloc(ptr->glb_selection_lns_gr5[ibin],tau_size*sizeof(double),ptr->error_message);
    class_alloc(ptr->glb_selection_lns_tau0_minus_tau_nc[ibin],tau_size*sizeof(double),ptr->error_message);

    //Get sampling values of tau_lensing
    class_call(transfer_selection_sampling_nc(ppr,pba,ppt,ptr,ibin,
					   tau0_minus_tau_lns_sources,tau_sources_size),
	       ptr->error_message,ptr->error_message);
    class_call(array_trapezoidal_mweights(tau0_minus_tau_lns_sources,tau_sources_size,
					  w_trapz_lns_sources,ptr->error_message),
	       ptr->error_message,ptr->error_message);
    class_call(transfer_selection_compute_nc(ppr,pba,ppt,ptr,selection,
					     tau0_minus_tau_lns_sources,
					     w_trapz_lns_sources,tau_sources_size,
					     pvecback,tau0,ibin),
	       ptr->error_message,ptr->error_message);

    //Get sampling values of tau
    class_call(transfer_lensing_sampling_nc(ppr,pba,ppt,ptr,ibin,tau0,
					 ptr->glb_selection_lns_tau0_minus_tau_nc[ibin],
					 ptr->glb_selection_lns_size_nc[ibin]),
	       ptr->error_message,ptr->error_message);

    int index_tau_sourcesb;
    class_alloc(array_bg,5*tau_sources_size*sizeof(double),ptr->error_message);
    
    for(index_tau_sourcesb=0;index_tau_sourcesb<tau_sources_size;index_tau_sourcesb++) {
      double z,e_bias,s_bias;

      class_call(background_at_tau(pba,tau0-tau0_minus_tau_lns_sources[index_tau_sourcesb],
				   pba->long_info,pba->inter_normal,&last_index,pvecback),
		 pba->error_message,ptr->error_message);
      array_bg[5*index_tau_sourcesb+1]=1./pvecback[pba->index_bg_a];
      array_bg[5*index_tau_sourcesb+2]=1./pvecback[pba->index_bg_H];
      array_bg[5*index_tau_sourcesb+3]=pvecback[pba->index_bg_H_prime];
      z=pba->a_today/pvecback[pba->index_bg_a]-1;
      
      if(ppt->has_nc_gr5) {
	class_call(transfer_interpolate_bias_function(ppr,pba,ppt,ptr,z,itr,
						      ptr->ez_z,ptr->ez_ez,ptr->ez_ddez,ptr->ez_size,
						      &last_index,&e_bias),
		   ptr->error_message,ptr->error_message);
      }
      else
	e_bias=0.;
      class_call(transfer_interpolate_bias_function(ppr,pba,ppt,ptr,z,itr,
						    ptr->sz_z,ptr->sz_sz,ptr->sz_ddsz,ptr->sz_size,
						    &last_index,&s_bias),
		 ptr->error_message,ptr->error_message);
      
      array_bg[5*index_tau_sourcesb+4]=e_bias;
      array_bg[5*index_tau_sourcesb+0]=s_bias;
    }
   
    for(index_tau=0;index_tau<tau_size;index_tau++) {
      double rescaling_lns=0,rescaling_gr4=0,rescaling_gr5=0;

      if(index_tau==tau_size-1) {
	if(ppt->has_nc_lens)
	  rescaling_lns=0;
	if(ppt->has_nc_gr4)
	  rescaling_gr4=0;
	if(ppt->has_nc_gr5)
	  rescaling_gr5=0;
      }
      else {
	int index_tau_sources;

	if(ppt->has_nc_lens)
	  rescaling_lns=0;
	if(ppt->has_nc_gr4)
	  rescaling_gr4=0;
	if(ppt->has_nc_gr5)
	  rescaling_gr5=0;
	
	for(index_tau_sources=0;index_tau_sources<tau_sources_size;index_tau_sources++) {
	  if((tau0_minus_tau_lns_sources[index_tau_sources]>0.) &&
	     (tau0_minus_tau_lns_sources[index_tau_sources]-
	      ptr->glb_selection_lns_tau0_minus_tau_nc[ibin][index_tau]>0.)) {
	    double s_bias=array_bg[5*index_tau_sources+0];

	    if(ppt->has_nc_lens) {
	      rescaling_lns+=(2.-5.*s_bias)/2.
		*(ptr->glb_selection_lns_tau0_minus_tau_nc[ibin][index_tau]-
		  tau0_minus_tau_lns_sources[index_tau_sources])
		/ptr->glb_selection_lns_tau0_minus_tau_nc[ibin][index_tau]
		/tau0_minus_tau_lns_sources[index_tau_sources]
		*selection[index_tau_sources]
		*w_trapz_lns_sources[index_tau_sources];
	    }

	    if(ppt->has_nc_gr4) {
	      rescaling_gr4+=(2.-5.*s_bias)
		/tau0_minus_tau_lns_sources[index_tau_sources]
		*selection[index_tau_sources]
		*w_trapz_lns_sources[index_tau_sources];
	    }

	    if(ppt->has_nc_gr5) {
	      double ia_h=array_bg[5*index_tau_sources+1];
	      double iH_h=array_bg[5*index_tau_sources+2];
	      double Hp_h=array_bg[5*index_tau_sources+3];
	      double e_bias=array_bg[5*index_tau_sources+4];

	      rescaling_gr5+=
		(1.+Hp_h*ia_h*iH_h*iH_h+(2.-5.*s_bias)*ia_h*iH_h/
		 tau0_minus_tau_lns_sources[index_tau_sources]+5.*s_bias-e_bias)
		*selection[index_tau_sources]
		*w_trapz_lns_sources[index_tau_sources];
	    }
	  }
	}
      }

      if(ppt->has_nc_lens)
	ptr->glb_selection_lns_lns[ibin][index_tau]=rescaling_lns;
      if(ppt->has_nc_gr4)
	ptr->glb_selection_lns_gr4[ibin][index_tau]=rescaling_gr4;
      if(ppt->has_nc_gr5)
	ptr->glb_selection_lns_gr5[ibin][index_tau]=rescaling_gr5;
    }
    free(array_bg);
  }

  free(pvecback);
  free(selection);
  free(tau0_minus_tau_lns_sources);
  free(w_trapz_lns_sources);

  return _SUCCESS_;
}

int init_selection_functions_nc_dam(
				 struct precision * ppr,
				 struct background * pba,
				 struct perturbs * ppt,
				 struct transfers * ptr
				 ) {
  int ibin;
  double tau0=pba->conformal_age;
  double *pvecback;

  class_alloc(pvecback,pba->bg_size*sizeof(double),ptr->error_message);
  class_alloc(ptr->glb_selection_0,ptr->selection_num_nc*sizeof(double *),ptr->error_message);
  class_alloc(ptr->glb_selection_0_tau0_minus_tau,ptr->selection_num_nc*sizeof(double *),ptr->error_message);
  class_alloc(ptr->glb_selection_0_size,ptr->selection_num_nc*sizeof(int),ptr->error_message);

  for(ibin=0;ibin<ptr->selection_num_nc;ibin++) {
    double tau_min,tau_mean,tau_max;
    double *w_trapz;
    int tau_size,l_limber;

    class_call(transfer_selection_times_nc(ppr,pba,ppt,ptr,ibin,&tau_min,&tau_mean,&tau_max),
	       ptr->error_message,ptr->error_message);

    //case selection=dirac;
    if(tau_min==tau_max) {
      tau_size=1;}
    else {
      //check that selection is well sampled
      tau_size=(int)ppr->selection_sampling;

      //value of l at which the code switches to limber approximation
      l_limber=ppr->l_switch_limber_for_cl_density;

      //check thab bessel is well sampled
      tau_size=MAX(tau_size,(int)((tau_max-tau_min)/((tau0-tau_mean)/l_limber))*ppr->selection_sampling_bessel);
    }

    ptr->glb_selection_0_size[ibin]=tau_size;

    class_alloc(ptr->glb_selection_0[ibin],tau_size*sizeof(double),ptr->error_message);
    class_alloc(ptr->glb_selection_0_tau0_minus_tau[ibin],tau_size*sizeof(double),ptr->error_message);
    class_alloc(w_trapz,tau_size*sizeof(double),ptr->error_message);

    class_call(transfer_selection_sampling_nc(ppr,pba,ppt,ptr,ibin,ptr->glb_selection_0_tau0_minus_tau[ibin],tau_size),
	       ptr->error_message,ptr->error_message);

    class_call(array_trapezoidal_mweights(ptr->glb_selection_0_tau0_minus_tau[ibin],tau_size,
					  w_trapz,ptr->error_message),ptr->error_message,ptr->error_message);

    class_call(transfer_selection_compute_nc(ppr,pba,ppt,ptr,
					     ptr->glb_selection_0[ibin],
					     ptr->glb_selection_0_tau0_minus_tau[ibin],
					     w_trapz,tau_size,pvecback,tau0,ibin),
	       ptr->error_message,ptr->error_message);

    free(w_trapz);
  }
  free(pvecback);

  return _SUCCESS_;
}

int init_selection_functions_wl_dam(
				    struct precision * ppr,
				    struct background * pba,
				    struct perturbs * ppt,
				    struct transfers * ptr
				    ) {
  int ibin;
  double tau0=pba->conformal_age;
  double *pvecback;

  class_alloc(pvecback,pba->bg_size*sizeof(double),ptr->error_message);
  class_alloc(ptr->glb_selection_0_wl,ptr->selection_num_wl*sizeof(double *),ptr->error_message);
  class_alloc(ptr->glb_selection_0_tau0_minus_tau_wl,ptr->selection_num_wl*sizeof(double *),ptr->error_message);
  class_alloc(ptr->glb_selection_0_size_wl,ptr->selection_num_wl*sizeof(int),ptr->error_message);

  for(ibin=0;ibin<ptr->selection_num_wl;ibin++) {
    double tau_min,tau_mean,tau_max;
    double *w_trapz;
    int tau_size,l_limber;

    class_call(transfer_selection_times_wl(ppr,pba,ppt,ptr,ibin,&tau_min,&tau_mean,&tau_max),
	       ptr->error_message,ptr->error_message);

    //case selection=dirac;
    if(tau_min==tau_max) {
      tau_size=1;}
    else {
      //check that selection is well sampled
      tau_size=(int)ppr->selection_sampling;

      //value of l at which the code switches to limber approximation
      l_limber=ppr->l_switch_limber_for_cl_lensing;

      //check thab bessel is well sampled
      tau_size=MAX(tau_size,(int)((tau_max-tau_min)/((tau0-tau_mean)/l_limber))*ppr->selection_sampling_bessel);
    }

    ptr->glb_selection_0_size_wl[ibin]=tau_size;

    class_alloc(ptr->glb_selection_0_wl[ibin],tau_size*sizeof(double),ptr->error_message);
    class_alloc(ptr->glb_selection_0_tau0_minus_tau_wl[ibin],tau_size*sizeof(double),ptr->error_message);
    class_alloc(w_trapz,tau_size*sizeof(double),ptr->error_message);

    class_call(transfer_selection_sampling_wl(ppr,pba,ppt,ptr,ibin,ptr->glb_selection_0_tau0_minus_tau_wl[ibin],tau_size),
	       ptr->error_message,ptr->error_message);

    class_call(array_trapezoidal_mweights(ptr->glb_selection_0_tau0_minus_tau_wl[ibin],tau_size,
					  w_trapz,ptr->error_message),ptr->error_message,ptr->error_message);

    class_call(transfer_selection_compute_wl(ppr,pba,ppt,ptr,
					  ptr->glb_selection_0_wl[ibin],
					  ptr->glb_selection_0_tau0_minus_tau_wl[ibin],
					  w_trapz,tau_size,pvecback,tau0,ibin),
	       ptr->error_message,ptr->error_message);

    free(w_trapz);
  }
  free(pvecback);

  return _SUCCESS_;
}

double get_f_nl_at_k(double k,struct perturbs * ppt)
{
  //Linear interpolation (shit?)
  int index_md=ppt->index_md_scalars;

  if(k<=ppt->k[index_md][0]) {
    double k0=ppt->k[index_md][0];
    return ppt->f_nl_kdep[0]*k0*k0/(k*k);
  }
  else {
    int index_k=0;
    double u;

    if(k>=ppt->k[index_md][ppt->k_size[index_md]-1])
      index_k=ppt->k_size[index_md]-1;
    else {
      index_k=0;
      while(ppt->k[index_md][index_k]<k)
	index_k++;
    }

    u=(k-ppt->k[index_md][index_k-1])/
      (ppt->k[index_md][index_k]-ppt->k[index_md][index_k-1]);

    return ppt->f_nl_kdep[index_k-1]*(1-u)+ppt->f_nl_kdep[index_k]*u;
  }
}

/**
 * Transfer function \f$ \Delta_l^{X} (q) \f$ at a given wavenumber q.
 *
 * For a given mode (scalar, vector, tensor), initial condition, type
 * (temperature, polarization, lensing, etc) and multipole, computes
 * the transfer function for an arbitary value of q by interpolating
 * between pre-computed values of q. This
 * function can be called from whatever module at whatever time,
 * provided that transfer_init() has been called before, and
 * transfer_free() has not been called yet.
 *
 * Wavenumbers are called q in this module and k in the perturbation
 * module. In flat universes k=q. In non-flat universes q and k differ
 * through q2 = k2 + K(1+m), where m=0,1,2 for scalar, vector,
 * tensor. q should be used throughout the transfer module, excepted
 * when interpolating or manipulating the source functions S(k,tau)
 * calculated in the perturbation module: for a given value of q, this
 * should be done at the corresponding k(q).
 *
 * @param index_md Input: index of requested mode
 * @param index_ic   Input: index of requested initial condition
 * @param index_tt   Input: index of requested type
 * @param index_l    Input: index of requested multipole
 * @param k          Input: any wavenumber
 * @param transfer_function Output: transfer function
 * @return the error status
 */

int transfer_functions_at_q(
                            struct transfers * ptr,
                            int index_md,
                            int index_ic,
                            int index_tt,
                            int index_l,
                            double q,
                            double * transfer_function
                            ) {
  /** Summary: */

  /** - interpolate in pre-computed table using array_interpolate_two() */
  class_call(array_interpolate_two(
                                   ptr->q,
                                   1,
                                   0,
                                   ptr->transfer[index_md]
                                   +((index_ic * ptr->tt_size[index_md] + index_tt) * ptr->l_size[index_md] + index_l)
                                   * ptr->q_size,
                                   1,
                                   ptr->q_size,
                                   q,
                                   transfer_function,
                                   1,
                                   ptr->error_message),
             ptr->error_message,
             ptr->error_message);

  return _SUCCESS_;
}

/**
 * This routine initializes the transfers structure, (in particular,
 * computes table of transfer functions \f$ \Delta_l^{X} (q) \f$)
 *
 * Main steps:
 *
 * - initialize all indices in the transfers structure
 *   and allocate all its arrays using transfer_indices_of_transfers().
 *
 * - for each thread (in case of parallel run), initialize the fields of a memory zone called the transfer_workspace with transfer_workspace_init()
 *
 * - loop over q values. For each q, compute the bessel functions if needed with transfer_update_HIS(), and defer the calculation of all transfer functions to transfer_compute_for_each_q()
 * - for each thread, free the the workspace with transfer_workspace_free()
 *
 * @param ppr Input : pointer to precision structure
 * @param pba Input : pointer to background structure
 * @param pth Input : pointer to thermodynamics structure
 * @param ppt Input : pointer to perturbation structure
 * @param ptr Output: pointer to initialized transfers structure
 * @return the error status
 */

int transfer_init(
                  struct precision * ppr,
                  struct background * pba,
                  struct thermo * pth,
                  struct perturbs * ppt,
                  struct nonlinear * pnl,
                  struct transfers * ptr
                  ) {

  /** Summary: */

  /** - define local variables */

  /* running index for wavenumbers */
  int index_q;

  /* conformal time today */
  double tau0;
  /* conformal time at recombination */
  double tau_rec;
  /* order of magnitude of the oscillation period of transfer functions */
  double q_period;

  /* maximum number of sampling times for transfer sources */
  int tau_size_max;

  /* array of sources S(k,tau), just taken from perturbation module,
     or transformed if non-linear corrections are needed
     sources[index_md][index_ic * ppt->tp_size[index_md] + index_tp][index_tau * ppt->k_size[index_md] + index_k]
  */
  double *** sources;

  /* array of source derivatives S''(k,tau)
     (second derivative with respect to k, not tau!),
     used to interpolate sources at the right values of k,
     sources_spline[index_md][index_ic * ppt->tp_size[index_md] + index_tp][index_tau * ppt->k_size[index_md] + index_k]
  */
  double *** sources_spline;

  /* pointer on workspace (one per thread if openmp) */
  struct transfer_workspace * ptw;

  /** - array with the correspondance between the index of sources in
      the perturbation module and in the transfer module,
      tp_of_tt[index_md][index_tt]
  */
  int ** tp_of_tt;

  /* structure containing the flat spherical bessel functions */

  HyperInterpStruct BIS;
  double xmax;

  /* This code can be optionally compiled with the openmp option for parallel computation.
     Inside parallel regions, the use of the command "return" is forbidden.
     For error management, instead of "return _FAILURE_", we will set the variable below
     to "abort = _TRUE_". This will lead to a "return _FAILURE_" jus after leaving the
     parallel region. */
  int abort;

#ifdef _OPENMP

  /* instrumentation times */
  double tstart, tstop, tspent;

#endif

  /** check whether any spectrum in harmonic space (i.e., any C_l's) is actually requested */

  if (ppt->has_cls == _FALSE_) {
    ptr->has_cls = _FALSE_;
    if (ptr->transfer_verbose > 0)
      printf("No harmonic space transfer functions to compute. Transfer module skipped.\n");
    return _SUCCESS_;
  }
  else
    ptr->has_cls = _TRUE_;

  if (ptr->transfer_verbose > 0)
    fprintf(stdout,"Computing transfers\n");

  /** get number of modes (scalars, tensors...) */

  ptr->md_size = ppt->md_size;

  /** - get conformal age / recombination time
      from background / thermodynamics structures
      (only place where these structures are used in this module) */

  tau0 = pba->conformal_age;
  tau_rec = pth->tau_rec;

  /** - correspondance between k and l depend on angular diameter
      diatance, i.e. on curvature. */

  ptr->angular_rescaling = pth->angular_rescaling;

  /** order of magnitude of the oscillation period of transfer functions */

  q_period = 2.*_PI_/(tau0-tau_rec)*ptr->angular_rescaling;

  /** - initialize all indices in the transfers structure and
      allocate all its arrays using transfer_indices_of_transfers() */

  class_call(transfer_indices_of_transfers(ppr,ppt,ptr,q_period,pba->K,pba->sgnK),
             ptr->error_message,
             ptr->error_message);

  /** - copy sources to a local array sources (in fact, only the pointers are copied, not the data), and eventually apply non-linear corrections to the sources */

  class_alloc(sources,
              ptr->md_size*sizeof(double**),
              ptr->error_message);

  class_call(transfer_perturbation_copy_sources_and_nl_corrections(ppt,pnl,ptr,sources),
             ptr->error_message,
             ptr->error_message);

  /** - spline all the sources passed by the perturbation module with respect to k (in order to interpolate later at a given value of k) */

  class_alloc(sources_spline,
              ptr->md_size*sizeof(double**),
              ptr->error_message);

  class_call(transfer_perturbation_source_spline(ppt,ptr,sources,sources_spline),
             ptr->error_message,
             ptr->error_message);

  /** - allocate and fill array describing the correspondence between perturbation types and transfer types */

  class_alloc(tp_of_tt,
              ptr->md_size*sizeof(int*),
              ptr->error_message);

  class_call(transfer_get_source_correspondence(ppt,ptr,tp_of_tt),
             ptr->error_message,
             ptr->error_message);

  /** - evaluate maximum number of sampled times in the transfer
      sources: needs to be known here, in order to allocate a large
      enough workspace */

  class_call(transfer_source_tau_size_max(ppr,pba,ppt,ptr,tau_rec,tau0,&tau_size_max),
             ptr->error_message,
             ptr->error_message);

  /** - compute flat spherical bessel functions */

  xmax = ptr->q[ptr->q_size-1]*tau0;
  if (pba->sgnK == -1)
    xmax *= (ptr->l[ptr->l_size_max-1]/ppr->hyper_flat_approximation_nu)/asinh(ptr->l[ptr->l_size_max-1]/ppr->hyper_flat_approximation_nu)*1.01;

  class_call(hyperspherical_HIS_create(0,
                                       1.,
                                       ptr->l_size_max,
                                       ptr->l,
                                       ppr->hyper_x_min,
                                       xmax,
                                       ppr->hyper_sampling_flat,
                                       ptr->l[ptr->l_size_max-1]+1,
                                       ppr->hyper_phi_min_abs,
                                       &BIS,
                                       ptr->error_message),
             ptr->error_message,
             ptr->error_message);

  /** - eventually read the selection and evolution functions */

  class_call(transfer_global_selection_read(ppt,ptr),
             ptr->error_message,
             ptr->error_message);

  //if there's nc, precompute selection functions
  if((ppt->has_nc_density) ||
     (ppt->has_nc_rsd1) ||
     (ppt->has_nc_rsd2) ||
     (ppt->has_nc_rsd3) ||
     (ppt->has_nc_gr1) ||
     (ppt->has_nc_gr2) ||
     (ppt->has_nc_gr3)) {
    class_call(init_selection_functions_nc_dam(ppr,pba,ppt,ptr),
	       ptr->error_message,ptr->error_message);
  }
  if((ppt->has_nc_lens) ||
     (ppt->has_nc_gr4) ||
     (ppt->has_nc_gr5)) {
    class_call(init_selection_functions_lensing_nc_dam(ppr,pba,ppt,ptr),
	       ptr->error_message,ptr->error_message);
  }
  if(ppt->has_weak_lensing) {
    if(ppt->has_lensing_shear) {
      class_call(init_selection_functions_lensing_wl_dam(ppr,pba,ppt,ptr),
		 ptr->error_message,ptr->error_message);
    }
    if(ppt->has_intrinsic_alignment) {
      class_call(init_selection_functions_wl_dam(ppr,pba,ppt,ptr),
		 ptr->error_message,ptr->error_message);
    }
  }
  /** (a.3.) workspace, allocated in a parallel zone since in openmp
      version there is one workspace per thread */

  /* initialize error management flag */
  abort = _FALSE_;

  /* beginning of parallel region */

#pragma omp parallel                                                    \
  shared(tau_size_max,ptr,ppr,pba,ppt,tp_of_tt,tau_rec,sources_spline,abort,BIS,tau0) \
  private(ptw,index_q,tstart,tstop,tspent)
  {

#ifdef _OPENMP
    tspent = 0.;
#endif

    /* allocate workspace */

    class_call_parallel(transfer_workspace_init(ptr,
                                                ppr,
                                                &ptw,
                                                ppt->tau_size,
                                                tau_size_max,
                                                pba->K,
                                                pba->sgnK,
                                                tau0-pth->tau_cut,
                                                &BIS),
                        ptr->error_message,
                        ptr->error_message);

    /** - loop over all wavenumbers (parallelised). For each wavenumber: */

#pragma omp for schedule (dynamic)

    for (index_q = 0; index_q < ptr->q_size; index_q++) {

#ifdef _OPENMP
      tstart = omp_get_wtime();
#endif

      if (ptr->transfer_verbose > 2)
        printf("Compute transfer for wavenumber [%d/%zu]\n",index_q,ptr->q_size-1);

      /* Update interpolation structure: */
      class_call_parallel(transfer_update_HIS(ppr,
                                              ptr,
                                              ptw,
                                              index_q,
                                              tau0),
                          ptr->error_message,
                          ptr->error_message);

      class_call_parallel(transfer_compute_for_each_q(ppr,
                                                      pba,
                                                      ppt,
                                                      ptr,
                                                      tp_of_tt,
                                                      index_q,
                                                      tau_size_max,
                                                      tau_rec,
                                                      sources,
                                                      sources_spline,
                                                      ptw),
                          ptr->error_message,
                          ptr->error_message);
      class_call_parallel(transfer_gather_for_each_q(ppr,
						     pba,
						     ppt,
						     ptr,
						     index_q),
                          ptr->error_message,
                          ptr->error_message);

#ifdef _OPENMP
      tstop = omp_get_wtime();

      tspent += tstop-tstart;
#endif

#pragma omp flush(abort)

    } /* end of loop over wavenumber */

    /* free workspace allocated inside parallel zone */
    class_call_parallel(transfer_workspace_free(ptr,ptw),
                        ptr->error_message,
                        ptr->error_message);

#ifdef _OPENMP
    if (ptr->transfer_verbose>1)
      printf("In %s: time spent in parallel region (loop over k's) = %e s for thread %d\n",
             __func__,tspent,omp_get_thread_num());
#endif

  } /* end of parallel region */

  if (abort == _TRUE_) return _FAILURE_;

  /* finally, free arrays allocated outside parallel zone */

  class_call(transfer_perturbation_sources_spline_free(ppt,ptr,sources_spline),
             ptr->error_message,
             ptr->error_message);

  class_call(transfer_perturbation_sources_free(ppt,pnl,ptr,sources),
             ptr->error_message,
             ptr->error_message);

  class_call(transfer_free_source_correspondence(ptr,tp_of_tt),
             ptr->error_message,
             ptr->error_message);

  class_call(hyperspherical_HIS_free(&BIS,ptr->error_message),
             ptr->error_message,
             ptr->error_message);

  if((ppt->has_cl_number_count) || (ppt->has_weak_lensing)) {
    free_selection_functions_dam(ppt,ptr);
  }

  return _SUCCESS_;
}

/**
 * This routine frees all the memory space allocated by transfer_init().
 *
 * To be called at the end of each run, only when no further calls to
 * transfer_functions_at_k() are needed.
 *
 * @param ptr Input: pointer to transfers structure (which fields must be freed)
 * @return the error status
 */

int transfer_free(
                  struct transfers * ptr
                  ) {

  int index_md,ii;

  if (ptr->has_cls == _TRUE_) {

    for (index_md = 0; index_md < ptr->md_size; index_md++) {
      free(ptr->l_size_tt[index_md]);
      free(ptr->l_size_tt_all[index_md]);
      free(ptr->transfer[index_md]);
      free(ptr->transfer_all[index_md]);
      free(ptr->k[index_md]);
    }

    free(ptr->tt_size);
    free(ptr->tt_size_all);
    free(ptr->l_size_tt);
    free(ptr->l_size_tt_all);
    free(ptr->l_size);
    free(ptr->l);
    free(ptr->q);
    free(ptr->k);
    free(ptr->transfer);
    free(ptr->transfer_all);

    for(ii=0;ii<ptr->n_tracers_wl;ii++) {
      if (ptr->nz_size_wl[ii] > 0) {
	free(ptr->nz_z_wl[ii]);
	free(ptr->nz_nz_wl[ii]);
	free(ptr->nz_ddnz_wl[ii]);
      }
      if (ptr->aI_size[ii] > 0) {
	free(ptr->aI_z[ii]);
	free(ptr->aI_aI[ii]);
	free(ptr->aI_ddaI[ii]);
      }
      if (ptr->fred_size[ii] > 0) {
	free(ptr->fred_z[ii]);
	free(ptr->fred_fred[ii]);
	free(ptr->fred_ddfred[ii]);
      }
    }

    for(ii=0;ii<ptr->n_tracers_nc;ii++) {
      if (ptr->nz_size_nc[ii] > 0) {
	free(ptr->nz_z_nc[ii]);
	free(ptr->nz_nz_nc[ii]);
	free(ptr->nz_ddnz_nc[ii]);
      }
      if (ptr->bz_size[ii] > 0) {
	free(ptr->bz_z[ii]);
	free(ptr->bz_bz[ii]);
	free(ptr->bz_ddbz[ii]);
      }
      if (ptr->sz_size[ii] > 0) {
	free(ptr->sz_z[ii]);
	free(ptr->sz_sz[ii]);
	free(ptr->sz_ddsz[ii]);
      }
      if (ptr->ez_size[ii] > 0) {
	free(ptr->ez_z[ii]);
	free(ptr->ez_ez[ii]);
	free(ptr->ez_ddez[ii]);
      }
    }
  }

  return _SUCCESS_;

}

/**
 * This routine defines all indices and allocates all tables
 * in the transfers structure
 *
 * Compute list of (k, l) values, allocate and fill corresponding
 * arrays in the transfers structure. Allocate the array of transfer
 * function tables.
 *
 * @param ppr Input : pointer to precision structure
 * @param ppt Input : pointer to perturbation structure
 * @param ptr Input/Output: pointer to transfer structure
 * @param rs_rec  Input : comoving distance to recombination
 * @return the error status
 */

int transfer_indices_of_transfers(
                                  struct precision * ppr,
                                  struct perturbs * ppt,
                                  struct transfers * ptr,
                                  double q_period,
                                  double K,
                                  int sgnK
                                  ) {

  /** Summary: */

  /** - define local variables */

  int index_md,index_tt,index_tt_common;

  /** define indices for transfer types */

  class_alloc(ptr->tt_size,ptr->md_size * sizeof(int),ptr->error_message);
  class_alloc(ptr->tt_size_all,ptr->md_size * sizeof(int),ptr->error_message);

  /** - type indices common to scalars and tensors */

  index_tt = 0;

  class_define_index(ptr->index_tt_t2,ppt->has_cl_cmb_temperature, index_tt,1);
  class_define_index(ptr->index_tt_e, ppt->has_cl_cmb_polarization,index_tt,1);

  index_tt_common=index_tt;

  /** - type indices for scalars */

  if (ppt->has_scalars == _TRUE_) {

    index_tt = index_tt_common;

    class_define_index(ptr->index_tt_t0,       ppt->has_cl_cmb_temperature,      index_tt,1);
    class_define_index(ptr->index_tt_t1,       ppt->has_cl_cmb_temperature,      index_tt,1);
    class_define_index(ptr->index_tt_lcmb,     ppt->has_cl_cmb_lensing_potential,index_tt,1);
    class_define_index(ptr->index_tt_density,  ppt->has_nc_density,              index_tt,ptr->selection_num_nc);
    class_define_index(ptr->index_tt_rsd,      ppt->has_nc_rsd1,                 index_tt,ptr->selection_num_nc);
    class_define_index(ptr->index_tt_d0,       ppt->has_nc_rsd2,                 index_tt,ptr->selection_num_nc);
    class_define_index(ptr->index_tt_d1,       ppt->has_nc_rsd3,                 index_tt,ptr->selection_num_nc);
    class_define_index(ptr->index_tt_nc_lens,  ppt->has_nc_lens,                 index_tt,ptr->selection_num_nc);
    class_define_index(ptr->index_tt_nc_g1,    ppt->has_nc_gr1,                  index_tt,ptr->selection_num_nc);
    class_define_index(ptr->index_tt_nc_g2,    ppt->has_nc_gr2,                  index_tt,ptr->selection_num_nc);
    class_define_index(ptr->index_tt_nc_g3,    ppt->has_nc_gr3,                  index_tt,ptr->selection_num_nc);
    class_define_index(ptr->index_tt_nc_g4,    ppt->has_nc_gr4,                  index_tt,ptr->selection_num_nc);
    class_define_index(ptr->index_tt_nc_g5,    ppt->has_nc_gr5,                  index_tt,ptr->selection_num_nc);
    class_define_index(ptr->index_tt_lensing,  ppt->has_lensing_shear,           index_tt,ptr->selection_num_wl);
    class_define_index(ptr->index_tt_alignment,ppt->has_intrinsic_alignment,     index_tt,ptr->selection_num_wl);

    ptr->tt_size[ppt->index_md_scalars]=index_tt;

  }

  /** - type indices for vectors */

  if (ppt->has_vectors == _TRUE_) {

    index_tt = index_tt_common;

    class_define_index(ptr->index_tt_t1,ppt->has_cl_cmb_temperature, index_tt,1);
    class_define_index(ptr->index_tt_b, ppt->has_cl_cmb_polarization,index_tt,1);

    ptr->tt_size[ppt->index_md_vectors]=index_tt;

  }

  /** - type indices for tensors */

  if (ppt->has_tensors == _TRUE_) {

    index_tt = index_tt_common;

    class_define_index(ptr->index_tt_b, ppt->has_cl_cmb_polarization,index_tt,1);

    ptr->tt_size[ppt->index_md_tensors]=index_tt;

  }

  //Indices for final transfer functions
  index_tt = 0;
  class_define_index(ptr->index_tt_tall,ppt->has_cl_cmb_temperature, index_tt,1);
  class_define_index(ptr->index_tt_eall,ppt->has_cl_cmb_polarization,index_tt,1);
  index_tt_common=index_tt;
  if (ppt->has_scalars == _TRUE_) {
    index_tt = index_tt_common;
    class_define_index(ptr->index_tt_pall ,ppt->has_cl_cmb_lensing_potential,index_tt,1);
    class_define_index(ptr->index_tt_ncall,ppt->has_cl_number_count         ,index_tt,ptr->selection_num_nc);
    class_define_index(ptr->index_tt_wlall,ppt->has_weak_lensing            ,index_tt,ptr->selection_num_wl);
    ptr->tt_size_all[ppt->index_md_scalars]=index_tt;
  }
  if (ppt->has_vectors == _TRUE_) {
    index_tt = index_tt_common;
    class_define_index(ptr->index_tt_ball,ppt->has_cl_cmb_polarization,index_tt,1);
    ptr->tt_size_all[ppt->index_md_vectors]=index_tt;
  }
  if (ppt->has_tensors == _TRUE_) {
    index_tt = index_tt_common;
    class_define_index(ptr->index_tt_ball,ppt->has_cl_cmb_polarization,index_tt,1);
    ptr->tt_size_all[ppt->index_md_tensors]=index_tt;
  }

  /** - allocate arrays of (k, l) values and transfer functions */

  /* number of l values for each mode and type,
     l_size_tt[index_md][index_tt], and maximized for each mode,
     l_size[index_md] */

  class_alloc(ptr->l_size,ptr->md_size * sizeof(int),ptr->error_message);

  class_alloc(ptr->l_size_tt,ptr->md_size * sizeof(int *),ptr->error_message);
  class_alloc(ptr->l_size_tt_all,ptr->md_size * sizeof(int *),ptr->error_message);

  for (index_md = 0; index_md < ptr->md_size; index_md++) {
    class_alloc(ptr->l_size_tt[index_md],ptr->tt_size[index_md] * sizeof(int),ptr->error_message);
    class_alloc(ptr->l_size_tt_all[index_md],ptr->tt_size_all[index_md] * sizeof(int),ptr->error_message);
  }

  /* array (of array) of transfer functions for each mode, transfer[index_md] */

  class_alloc(ptr->transfer,ptr->md_size * sizeof(double *),ptr->error_message);
  class_alloc(ptr->transfer_all,ptr->md_size * sizeof(double *),ptr->error_message);

  /** get q values using transfer_get_q_list() */

  class_call(transfer_get_q_list(ppr,ppt,ptr,q_period,K,sgnK),
             ptr->error_message,
             ptr->error_message);

  /** get k values using transfer_get_k_list() */
  class_call(transfer_get_k_list(ppt,ptr,K),
             ptr->error_message,
             ptr->error_message);

  /** get l values using transfer_get_l_list() */
  class_call(transfer_get_l_list(ppr,ppt,ptr),
             ptr->error_message,
             ptr->error_message);

  /** - loop over modes (scalar, etc). For each mode: */

  for (index_md = 0; index_md < ptr->md_size; index_md++) {

    /** allocate arrays of transfer functions, (ptr->transfer[index_md])[index_ic][index_tt][index_l][index_k] */
    class_alloc(ptr->transfer[index_md],
                ppt->ic_size[index_md]*ptr->tt_size[index_md]*ptr->l_size[index_md]*ptr->q_size*sizeof(double),
                ptr->error_message);
#ifdef _CLASST_DEBUG
    unsigned long long d_tsize=sizeof(double)*ptr->tt_size[index_md]*
      ((unsigned long long)(ptr->l_size[index_md]*ptr->q_size));
    double d_tsize_GB=(double)d_tsize/(1024.*1024.*1024.);
    printf("Mode %d : %.3lf GB for transfer functions here\n",index_md,d_tsize_GB);
    printf("          tt_size = %d, l_size= %d, q_size= %d\n",
	   ptr->tt_size[index_md],ptr->l_size[index_md],(int)(ptr->q_size));
#endif //_CLASST_DEBUG
    class_alloc(ptr->transfer_all[index_md],
                ppt->ic_size[index_md] * ptr->tt_size_all[index_md] * ptr->l_size[index_md] *
		ptr->q_size * sizeof(double),
                ptr->error_message);
#ifdef _CLASST_DEBUG
    d_tsize=sizeof(double)*ptr->tt_size_all[index_md]*
      ((unsigned long long)(ptr->l_size[index_md]*ptr->q_size));
    d_tsize_GB=(double)d_tsize/(1024.*1024.*1024.);
    printf("Final transfers\n");
    printf("          %.3lf GB for transfer functions here\n",index_md,d_tsize_GB);
    printf("          tt_size = %d, l_size= %d, q_size= %d\n",
	   ptr->tt_size_all[index_md],ptr->l_size[index_md],(int)(ptr->q_size));
#endif //_CLASST_DEBUG
  }

  return _SUCCESS_;

}

int transfer_perturbation_copy_sources_and_nl_corrections(
                                                          struct perturbs * ppt,
                                                          struct nonlinear * pnl,
                                                          struct transfers * ptr,
                                                          double *** sources
                                                          ) {
  int index_md;
  int index_ic;
  int index_tp;
  int index_k;
  int index_tau;

  for (index_md = 0; index_md < ptr->md_size; index_md++) {

    class_alloc(sources[index_md],
                ppt->ic_size[index_md]*ppt->tp_size[index_md]*sizeof(double*),
                ptr->error_message);

    for (index_ic = 0; index_ic < ppt->ic_size[index_md]; index_ic++) {

      for (index_tp = 0; index_tp < ppt->tp_size[index_md]; index_tp++) {

        if ((pnl->method != nl_none) && (_scalars_) &&
            (((ppt->has_source_delta_m == _TRUE_) && (index_tp == ppt->index_tp_delta_m)) ||
             ((ppt->has_source_theta_m == _TRUE_) && (index_tp == ppt->index_tp_theta_m)) ||
             ((ppt->has_source_phi == _TRUE_) && (index_tp == ppt->index_tp_phi)) ||
             ((ppt->has_source_phi_prime == _TRUE_) && (index_tp == ppt->index_tp_phi_prime)) ||
             ((ppt->has_source_phi_plus_psi == _TRUE_) && (index_tp == ppt->index_tp_phi_plus_psi)) ||
             ((ppt->has_source_psi == _TRUE_) && (index_tp == ppt->index_tp_psi)))) {

          class_alloc(sources[index_md][index_ic * ppt->tp_size[index_md] + index_tp],
                      ppt->k_size[index_md]*ppt->tau_size*sizeof(double),
                      ptr->error_message);

          for (index_tau=0; index_tau<ppt->tau_size; index_tau++) {
            for (index_k=0; index_k<ppt->k_size[index_md]; index_k++) {
              sources[index_md]
                [index_ic * ppt->tp_size[index_md] + index_tp]
                [index_tau * ppt->k_size[index_md] + index_k] =
                ppt->sources[index_md]
                [index_ic * ppt->tp_size[index_md] + index_tp]
                [index_tau * ppt->k_size[index_md] + index_k]
                * pnl->nl_corr_density[index_tau * ppt->k_size[index_md] + index_k];
            }
          }
        }
        else {
          sources[index_md][index_ic * ppt->tp_size[index_md] + index_tp] =
            ppt->sources[index_md][index_ic * ppt->tp_size[index_md] + index_tp];
        }
      }
    }
  }

  return _SUCCESS_;

}


int transfer_perturbation_source_spline(
                                        struct perturbs * ppt,
                                        struct transfers * ptr,
                                        double *** sources,
                                        double *** sources_spline
                                        ) {
  int index_md;
  int index_ic;
  int index_tp;

  for (index_md = 0; index_md < ptr->md_size; index_md++) {

    class_alloc(sources_spline[index_md],
                ppt->ic_size[index_md]*ppt->tp_size[index_md]*sizeof(double*),
                ptr->error_message);

    for (index_ic = 0; index_ic < ppt->ic_size[index_md]; index_ic++) {

      for (index_tp = 0; index_tp < ppt->tp_size[index_md]; index_tp++) {

        class_alloc(sources_spline[index_md][index_ic * ppt->tp_size[index_md] + index_tp],
                    ppt->k_size[index_md]*ppt->tau_size*sizeof(double),
                    ptr->error_message);

        class_call(array_spline_table_columns2(ppt->k[index_md],
                                               ppt->k_size[index_md],
                                               sources[index_md][index_ic * ppt->tp_size[index_md] + index_tp],
                                               ppt->tau_size,
                                               sources_spline[index_md][index_ic * ppt->tp_size[index_md] + index_tp],
                                               _SPLINE_EST_DERIV_,
                                               ptr->error_message),
                   ptr->error_message,
                   ptr->error_message);

      }
    }
  }

  return _SUCCESS_;

}

int transfer_perturbation_sources_free(
                                       struct perturbs * ppt,
                                       struct nonlinear * pnl,
                                       struct transfers * ptr,
                                       double *** sources
                                       ) {
  int index_md;
  int index_ic;
  int index_tp;

  for (index_md = 0; index_md < ptr->md_size; index_md++) {
    for (index_ic = 0; index_ic < ppt->ic_size[index_md]; index_ic++) {
      for (index_tp = 0; index_tp < ppt->tp_size[index_md]; index_tp++) {
        if ((pnl->method != nl_none) && (_scalars_) &&
            (((ppt->has_source_delta_m == _TRUE_) && (index_tp == ppt->index_tp_delta_m)) ||
             ((ppt->has_source_theta_m == _TRUE_) && (index_tp == ppt->index_tp_theta_m)) ||
             ((ppt->has_source_phi == _TRUE_) && (index_tp == ppt->index_tp_phi)) ||
             ((ppt->has_source_phi_prime == _TRUE_) && (index_tp == ppt->index_tp_phi_prime)) ||
             ((ppt->has_source_phi_plus_psi == _TRUE_) && (index_tp == ppt->index_tp_phi_plus_psi)) ||
             ((ppt->has_source_psi == _TRUE_) && (index_tp == ppt->index_tp_psi)))) {

          free(sources[index_md][index_ic * ppt->tp_size[index_md] + index_tp]);
        }
      }
    }
    free(sources[index_md]);
  }
  free(sources);

  return _SUCCESS_;
}

int transfer_perturbation_sources_spline_free(
                                              struct perturbs * ppt,
                                              struct transfers * ptr,
                                              double *** sources_spline
                                              ) {
  int index_md;
  int index_ic;
  int index_tp;

  for (index_md = 0; index_md < ptr->md_size; index_md++) {
    for (index_ic = 0; index_ic < ppt->ic_size[index_md]; index_ic++) {
      for (index_tp = 0; index_tp < ppt->tp_size[index_md]; index_tp++) {
        free(sources_spline[index_md][index_ic * ppt->tp_size[index_md] + index_tp]);
      }
    }
    free(sources_spline[index_md]);
  }
  free(sources_spline);

  return _SUCCESS_;
}

/**
 * This routine defines the number and values of mutipoles l for all modes.
 *
 * @param ppr  Input : pointer to precision structure
 * @param ppt  Input : pointer to perturbation structure
 * @param ptr  Input/Output : pointer to transfers structure containing l's
 * @return the error status
 */

int transfer_get_l_list(
                        struct precision * ppr,
                        struct perturbs * ppt,
                        struct transfers * ptr
                        ) {

  int index_l;
  int l_max=0;
  int index_md;
  int index_tt;
  int increment,current_l;

  /* check that largests need value of l_max */

  if (ppt->has_cls == _TRUE_) {

    if (ppt->has_scalars == _TRUE_) {

      if ((ppt->has_cl_cmb_temperature == _TRUE_) ||
          (ppt->has_cl_cmb_polarization == _TRUE_) ||
          (ppt->has_cl_cmb_lensing_potential == _TRUE_))
        l_max=MAX(ppt->l_scalar_max,l_max);

      if ((ppt->has_weak_lensing == _TRUE_) ||
          (ppt->has_cl_number_count == _TRUE_))
        l_max=MAX(ppt->l_lss_max,l_max);
    }

    if (ppt->has_tensors == _TRUE_)
      l_max=MAX(ppt->l_tensor_max,l_max);

  }

  /* allocate and fill l array */

  /** - start from l = 2 and increase with logarithmic step */

  index_l = 0;
  current_l = 2;
  increment = MAX((int)(current_l * (pow(ppr->l_logstep,ptr->angular_rescaling)-1.)),1);

  while (((current_l+increment) < l_max) &&
         (increment < ppr->l_linstep*ptr->angular_rescaling)) {

    index_l ++;
    current_l += increment;
    increment = MAX((int)(current_l * (pow(ppr->l_logstep,ptr->angular_rescaling)-1.)),1);

  }

  /** - when the logarithmic step becomes larger than some linear step,
      stick to this linear step till l_max */

  increment = ppr->l_linstep*ptr->angular_rescaling;

  while ((current_l+increment) <= l_max) {

    index_l ++;
    current_l += increment;

  }

  /** - last value set to exactly l_max */

  if (current_l != l_max) {

    index_l ++;
    current_l = l_max;

  }

  ptr->l_size_max = index_l+1;

  /** - so far we just counted the number of values. Now repeat the
      whole thing but fill array with values. */

  class_alloc(ptr->l,ptr->l_size_max*sizeof(int),ptr->error_message);

  index_l = 0;
  ptr->l[0] = 2;
  increment = MAX((int)(ptr->l[0] * (pow(ppr->l_logstep,ptr->angular_rescaling)-1.)),1);

  while (((ptr->l[index_l]+increment) < l_max) &&
         (increment < ppr->l_linstep*ptr->angular_rescaling)) {

    index_l ++;
    ptr->l[index_l]=ptr->l[index_l-1]+increment;
    increment = MAX((int)(ptr->l[index_l] * (pow(ppr->l_logstep,ptr->angular_rescaling)-1.)),1);

  }

  increment = ppr->l_linstep*ptr->angular_rescaling;

  while ((ptr->l[index_l]+increment) <= l_max) {

    index_l ++;
    ptr->l[index_l]=ptr->l[index_l-1]+increment;

  }

  if (ptr->l[index_l] != l_max) {

    index_l ++;
    ptr->l[index_l]= l_max;

  }

  /* for each mode and type, find relevant size of l array,
     l_size_tt[index_md][index_tt] (since for some modes and types
     l_max can be smaller). Also, maximize this size for each mode to
     find l_size[index_md]. */

  for (index_md=0; index_md < ppt->md_size; index_md++) {

    ptr->l_size[index_md] = 0;

    for (index_tt=0;index_tt<ptr->tt_size[index_md];index_tt++) {

      if (_scalars_) {
        if ((ppt->has_cl_cmb_temperature==_TRUE_) &&
            ((index_tt==ptr->index_tt_t0)||(index_tt==ptr->index_tt_t1)||(index_tt==ptr->index_tt_t2)))
          l_max=ppt->l_scalar_max;

        if ((ppt->has_cl_cmb_polarization == _TRUE_) && (index_tt == ptr->index_tt_e))
          l_max=ppt->l_scalar_max;

        if ((ppt->has_cl_cmb_lensing_potential == _TRUE_) && (index_tt == ptr->index_tt_lcmb))
          l_max=ppt->l_scalar_max;

        if ((_index_tt_in_range_(ptr->index_tt_density, ptr->selection_num_nc, ppt->has_nc_density)) ||
            (_index_tt_in_range_(ptr->index_tt_rsd,     ptr->selection_num_nc, ppt->has_nc_rsd1)) ||
            (_index_tt_in_range_(ptr->index_tt_d0,      ptr->selection_num_nc, ppt->has_nc_rsd2)) ||
            (_index_tt_in_range_(ptr->index_tt_d1,      ptr->selection_num_nc, ppt->has_nc_rsd3)) ||
            (_index_tt_in_range_(ptr->index_tt_nc_lens, ptr->selection_num_nc, ppt->has_nc_lens))||
            (_index_tt_in_range_(ptr->index_tt_nc_g1,   ptr->selection_num_nc, ppt->has_nc_gr1))  ||
            (_index_tt_in_range_(ptr->index_tt_nc_g2,   ptr->selection_num_nc, ppt->has_nc_gr2))  ||
            (_index_tt_in_range_(ptr->index_tt_nc_g3,   ptr->selection_num_nc, ppt->has_nc_gr3))  ||
            (_index_tt_in_range_(ptr->index_tt_nc_g4,   ptr->selection_num_nc, ppt->has_nc_gr4))  ||
            (_index_tt_in_range_(ptr->index_tt_nc_g5,   ptr->selection_num_nc, ppt->has_nc_gr5))
            )
          l_max=ppt->l_lss_max;

        if ((_index_tt_in_range_(ptr->index_tt_lensing,  ptr->selection_num_wl,ppt->has_lensing_shear)) ||
	    (_index_tt_in_range_(ptr->index_tt_alignment,ptr->selection_num_wl,ppt->has_intrinsic_alignment))
	    )
          l_max=ppt->l_lss_max;

      }

      if (_tensors_) {
        l_max = ppt->l_tensor_max;
      }

      class_test(l_max > ptr->l[ptr->l_size_max-1],
                 ptr->error_message,
                 "For mode %d, type %d, asked for l_max=%d greater than in Bessel table where l_max=%d",
                 index_md,
                 index_tt,
                 l_max,
                 ptr->l[ptr->l_size_max-1]);

      index_l=0;
      while (ptr->l[index_l] < l_max) index_l++;
      ptr->l_size_tt[index_md][index_tt]=index_l+1;

      if (ptr->l_size_tt[index_md][index_tt] < ptr->l_size_max)
        ptr->l_size_tt[index_md][index_tt]++;
      if (ptr->l_size_tt[index_md][index_tt] < ptr->l_size_max)
        ptr->l_size_tt[index_md][index_tt]++;

      ptr->l_size[index_md] = MAX(ptr->l_size[index_md],ptr->l_size_tt[index_md][index_tt]);

    }

    for (index_tt=0;index_tt<ptr->tt_size_all[index_md];index_tt++) {
      if (_scalars_) {
        if ((ppt->has_cl_cmb_temperature==_TRUE_) && (index_tt==ptr->index_tt_tall))
          l_max=ppt->l_scalar_max;

        if ((ppt->has_cl_cmb_polarization == _TRUE_) && (index_tt == ptr->index_tt_eall))
          l_max=ppt->l_scalar_max;

        if ((ppt->has_cl_cmb_lensing_potential == _TRUE_) && (index_tt == ptr->index_tt_pall))
          l_max=ppt->l_scalar_max;

        if ((_index_tt_in_range_(ptr->index_tt_ncall,ptr->selection_num_nc,ppt->has_cl_number_count)) ||
            (_index_tt_in_range_(ptr->index_tt_wlall,ptr->selection_num_wl,ppt->has_weak_lensing)))
          l_max=ppt->l_lss_max;
      }
      if (_tensors_) {
        l_max = ppt->l_tensor_max;
      }

      class_test(l_max > ptr->l[ptr->l_size_max-1],
                 ptr->error_message,
                 "For mode %d, type %d, asked for l_max=%d greater than in Bessel table where l_max=%d",
                 index_md,
                 index_tt,
                 l_max,
                 ptr->l[ptr->l_size_max-1]);

      index_l=0;
      while (ptr->l[index_l] < l_max) index_l++;
      ptr->l_size_tt_all[index_md][index_tt]=index_l+1;

      if (ptr->l_size_tt_all[index_md][index_tt] < ptr->l_size_max)
        ptr->l_size_tt_all[index_md][index_tt]++;
      if (ptr->l_size_tt_all[index_md][index_tt] < ptr->l_size_max)
        ptr->l_size_tt_all[index_md][index_tt]++;
    }
  }

  return _SUCCESS_;

}

/**
 * This routine defines the number and values of wavenumbers q for
 * each mode (goes smoothly from logarithmic step for small q's to
 * linear step for large q's).
 *
 * @param ppr     Input : pointer to precision structure
 * @param ppt     Input : pointer to perturbation structure
 * @param ptr     Input/Output : pointer to transfers structure containing q's
 * @param rs_rec  Input : comoving distance to recombination
 * @param index_md Input: index of requested mode (scalar, tensor, etc)
 * @return the error status
 */

int transfer_get_q_list(
                        struct precision * ppr,
                        struct perturbs * ppt,
                        struct transfers * ptr,
                        double q_period,
                        double K,
                        int sgnK
                        ) {

  int index_q;
  double q,q_min=0.,q_max=0.,q_step,k_max;
  int nu, nu_min, nu_proposed;
  int q_size_max;
  double q_approximation;
  double last_step=0.;
  int last_index=0;
  double q_logstep_spline;
  double q_logstep_trapzd;
  int index_md;

  /* first and last value in flat case*/

  if (sgnK == 0) {
    q_min = ppt->k_min;

    q_max = 0.;
    for (index_md=0; index_md<ppt->md_size; index_md++) {
      q_max = MAX(q_max,ppt->k[index_md][ppt->k_size_cl[index_md]-1]);
    }

    K=0;
  }

  /* first and last value in open case*/

  else if (sgnK == -1) {
    q_min = sqrt(ppt->k_min*ppt->k_min+K);

    k_max = 0.;
    for (index_md=0; index_md<ppt->md_size; index_md++) {
      k_max = MAX(k_max,ppt->k[index_md][ppt->k_size_cl[index_md]-1]);
    }

    q_max = sqrt(k_max*k_max+K);
    if (ppt->has_vectors == _TRUE_)
      q_max = MIN(q_max,sqrt(k_max*k_max+2.*K));
    if (ppt->has_tensors == _TRUE_)
      q_max = MIN(q_max,sqrt(k_max*k_max+3.*K));
  }

  /* first and last value in closed case*/

  else if (sgnK == 1) {
    nu_min = 3;
    q_min = nu_min * sqrt(K);

    q_max = 0.;
    for (index_md=0; index_md<ppt->md_size; index_md++) {
      q_max = MAX(q_max,ppt->k[index_md][ppt->k_size_cl[index_md]-1]);
    }
  }

  /* adjust the parameter governing the log step size to curvature */

  q_logstep_spline = ppr->q_logstep_spline/pow(ptr->angular_rescaling,ppr->q_logstep_open);
  q_logstep_trapzd = ppr->q_logstep_trapzd;

  /* very conservative estimate of number of values */

  if (sgnK == 1) {

    q_approximation = MIN(ppr->hyper_flat_approximation_nu,(q_max/sqrt(K)));

    /* max contribution from integer nu values */
    q_step = 1.+q_period*ppr->q_logstep_trapzd;
    q_size_max = 2*(int)(log(q_approximation/q_min)/log(q_step));

    q_step = q_period*ppr->q_linstep;
    q_size_max += 2*(int)((q_approximation-q_min)/q_step);

    /* max contribution from non-integer nu values */
    q_step = 1.+q_period*ppr->q_logstep_spline;
    q_size_max += 2*(int)(log(q_max/q_approximation)/log(q_step));

    q_step = q_period*ppr->q_linstep;
    q_size_max += 2*(int)((q_max-q_approximation)/q_step);

  }
  else {

    /* max contribution from non-integer nu values */
    q_step = 1.+q_period*ppr->q_logstep_spline;
    q_size_max = 5*(int)(log(q_max/q_min)/log(q_step));

    q_step = q_period*ppr->q_linstep;
    q_size_max += 5*(int)((q_max-q_min)/q_step);

  }

  /* create array with this conservative size estimate. The exact size
     will be readjusted below, after filling the array. */

  class_alloc(ptr->q,
              q_size_max*sizeof(double),
              ptr->error_message);

  /* assign the first value before starting the loop */

  index_q = 0;
  ptr->q[index_q] = q_min;
  nu = 3;
  index_q++;

  /* loop over the values */

  while (ptr->q[index_q-1] < q_max) {

    class_test(index_q >= q_size_max,ptr->error_message,"buggy q-list definition");

    /* step size formula in flat/open case. Step goes gradually from
       logarithmic to linear:

       - in the small q limit, it is logarithmic with: (delta q / q) =
       q_period * q_logstep_spline

       - in the large q limit, it is linear with: (delta q) = q_period
       * ppr->q_linstep
       */

    if (sgnK<=0) {

      q = ptr->q[index_q-1]
        + q_period * ppr->q_linstep * ptr->q[index_q-1]
        / (ptr->q[index_q-1] + ppr->q_linstep/q_logstep_spline);

    }

    /* step size formula in closed case. Same thing excepted that:

       - in the small q limit, the logarithmic step is reduced, being
       given by q_logstep_trapzd, and values are rounded to integer
       values of nu=q/sqrt(K). This happens as long as
       nu<nu_flat_approximation

       - for nu>nu_flat_approximation, the step gradually catches up
       the same expression as in the flat/opne case, and there is no
       need to round up to integer nu's.
    */

    else {

      if (nu < (int)ppr->hyper_flat_approximation_nu) {

        q = ptr->q[index_q-1]
          + q_period * ppr->q_linstep * ptr->q[index_q-1]
          / (ptr->q[index_q-1] + ppr->q_linstep/q_logstep_trapzd);

        nu_proposed = (int)(q/sqrt(K));
        if (nu_proposed <= nu+1)
          nu = nu+1;
        else
          nu = nu_proposed;

        q = nu*sqrt(K);
        last_step = q - ptr->q[index_q-1];
        last_index = index_q+1;
      }
      else {

        q_step = q_period * ppr->q_linstep * ptr->q[index_q-1] / (ptr->q[index_q-1] + ppr->q_linstep/q_logstep_spline);

        if (index_q-last_index < (int)ppr->q_numstep_transition)
          q = ptr->q[index_q-1] + (1-(double)(index_q-last_index)/ppr->q_numstep_transition) * last_step + (double)(index_q-last_index)/ppr->q_numstep_transition * q_step;
        else
          q = ptr->q[index_q-1] + q_step;
      }
    }

    ptr->q[index_q] = q;
    index_q++;
  }

  /* infer total number of values (also checking if we overshooted the last point) */

  if (ptr->q[index_q-1] > q_max)
    ptr->q_size=index_q-1;
  else
    ptr->q_size=index_q;

  class_test(ptr->q_size<2,ptr->error_message,"buggy q-list definition");

  /* now, readjust array size */

  class_realloc(ptr->q,
                ptr->q,
                ptr->q_size*sizeof(double),
                ptr->error_message);

  /* in curved universe, check at which index the flat rescaling
     approximation will start being used */

  if (sgnK != 0) {

    q_approximation = ppr->hyper_flat_approximation_nu * sqrt(sgnK*K);
    for (ptr->index_q_flat_approximation=0;
         ptr->index_q_flat_approximation < ptr->q_size-1;
         ptr->index_q_flat_approximation++) {
      if (ptr->q[ptr->index_q_flat_approximation] > q_approximation) break;
    }
    if (ptr->transfer_verbose > 1)
      printf("Flat bessel approximation spares hyperspherical bessel computations for %zu wavenumebrs over a total of %zu\n",
	     ptr->q_size-ptr->index_q_flat_approximation,ptr->q_size);
  }

  return _SUCCESS_;

}

/**
 * This routine infers from the q values a list of corresponding k
 * avlues for each mode.
 *
 * @param ppt     Input : pointer to perturbation structure
 * @param ptr     Input/Output : pointer to transfers structure containing q's
 * @param K       Input : spatial curvature
 * @return the error status
 */

int transfer_get_k_list(
                        struct perturbs * ppt,
                        struct transfers * ptr,
                        double K
                        ) {

  int index_md;
  int index_q;
  double m=0.;

  class_alloc(ptr->k,ptr->md_size*sizeof(double*),ptr->error_message);

  for (index_md = 0; index_md <  ptr->md_size; index_md++) {

    class_alloc(ptr->k[index_md],ptr->q_size*sizeof(double),ptr->error_message);

    if (_scalars_) {
      m=0.;
    }
    if (_vectors_) {
      m=1.;
    }
    if (_tensors_) {
      m=2.;
    }

    for (index_q=0; index_q < ptr->q_size; index_q++) {
      ptr->k[index_md][index_q] = sqrt(ptr->q[index_q]*ptr->q[index_q]-K*(m+1.));
    }

    if (ptr->k[index_md][0] < ppt->k[index_md][0]){
      /** If ptr->k[index_md][0] < ppt->k[index_md][0] at the level of rounding,
          adjust first value of k_list to avoid interpolation errors: */
      if ((ppt->k[index_md][0]-ptr->k[index_md][0]) < 10.*DBL_EPSILON){
        ptr->k[index_md][0] = ppt->k[index_md][0];
      }
      else{
        class_stop(ptr->error_message,
                   "bug in k_list calculation: in perturbation module k_min=%e, in transfer module k_min[mode=%d]=%e, interpolation impossible",
                   ppt->k[0][0],
                   index_md,
                   ptr->k[index_md][0]);
      }
    }

    /**
       class_test(ptr->k[index_md][0] < ppt->k[index_md][0],
       ptr->error_message,
       "bug in k_list calculation: in perturbation module k_min=%e, in transfer module k_min[mode=%d]=%e, interpolation impossible",
       ppt->k[0][0],
       index_md,
       ptr->k[index_md][0]);
    */
    class_test(ptr->k[index_md][ptr->q_size-1] > ppt->k[0][ppt->k_size_cl[0]-1],
               ptr->error_message,
               "bug in k_list calculation: in perturbation module k_max=%e, in transfer module k_max[mode=%d]=%e, interpolation impossible",
               ppt->k[0][ppt->k_size_cl[0]],
               index_md,
               ptr->k[index_md][ptr->q_size-1]);
  }

  return _SUCCESS_;

}

/**
 * This routine defines the correspondence between the sources in the
 * perturbation and transfer module.
 *
 * @param ppt  Input : pointer to perturbation structure
 * @param ptr  Input : pointer to transfers structure containing l's
 * @param index_md : Input: index of mode (scalar, tensor...)
 * @param tp_of_tt : Input/Output: array with the correspondance (allocated before, filled here)
 * @return the error status
 */

int transfer_get_source_correspondence(
                                       struct perturbs * ppt,
                                       struct transfers * ptr,
                                       int ** tp_of_tt
                                       ) {

  /* running index on modes */
  int index_md;

  /* running index on transfer types */
  int index_tt;

  /** - which source are we considering? Define correspondence
      between transfer types and source types */

  for (index_md = 0; index_md < ptr->md_size; index_md++) {

    class_alloc(tp_of_tt[index_md],ptr->tt_size[index_md]*sizeof(int),ptr->error_message);

    for (index_tt=0; index_tt<ptr->tt_size[index_md]; index_tt++) {

      if (_scalars_) {

        if ((ppt->has_cl_cmb_temperature == _TRUE_) && (index_tt == ptr->index_tt_t0))
          tp_of_tt[index_md][index_tt]=ppt->index_tp_t0;
        if ((ppt->has_cl_cmb_temperature == _TRUE_) && (index_tt == ptr->index_tt_t1))
          tp_of_tt[index_md][index_tt]=ppt->index_tp_t1;
        if ((ppt->has_cl_cmb_temperature == _TRUE_) && (index_tt == ptr->index_tt_t2))
          tp_of_tt[index_md][index_tt]=ppt->index_tp_t2;

        if ((ppt->has_cl_cmb_polarization == _TRUE_) && (index_tt == ptr->index_tt_e))
          tp_of_tt[index_md][index_tt]=ppt->index_tp_p;

        if ((ppt->has_cl_cmb_lensing_potential == _TRUE_) && (index_tt == ptr->index_tt_lcmb))
          tp_of_tt[index_md][index_tt]=ppt->index_tp_phi_plus_psi;

        if (_index_tt_in_range_(ptr->index_tt_density, ptr->selection_num_nc, ppt->has_nc_density))
          tp_of_tt[index_md][index_tt]=ppt->index_tp_delta_m;
        if (_index_tt_in_range_(ptr->index_tt_rsd,     ptr->selection_num_nc, ppt->has_nc_rsd1))
          tp_of_tt[index_md][index_tt]=ppt->index_tp_theta_m;
        if (_index_tt_in_range_(ptr->index_tt_d0,      ptr->selection_num_nc, ppt->has_nc_rsd2))
          tp_of_tt[index_md][index_tt]=ppt->index_tp_theta_m;
        if (_index_tt_in_range_(ptr->index_tt_d1,      ptr->selection_num_nc, ppt->has_nc_rsd3))
          tp_of_tt[index_md][index_tt]=ppt->index_tp_theta_m;
        if (_index_tt_in_range_(ptr->index_tt_nc_lens, ptr->selection_num_nc, ppt->has_nc_lens))
          tp_of_tt[index_md][index_tt]=ppt->index_tp_phi_plus_psi;
        if (_index_tt_in_range_(ptr->index_tt_nc_g1,   ptr->selection_num_nc, ppt->has_nc_gr1))
          tp_of_tt[index_md][index_tt]=ppt->index_tp_psi;
        if (_index_tt_in_range_(ptr->index_tt_nc_g2,   ptr->selection_num_nc, ppt->has_nc_gr2))
          tp_of_tt[index_md][index_tt]=ppt->index_tp_phi;
        if (_index_tt_in_range_(ptr->index_tt_nc_g3,   ptr->selection_num_nc, ppt->has_nc_gr3))
          tp_of_tt[index_md][index_tt]=ppt->index_tp_phi_prime;
        if (_index_tt_in_range_(ptr->index_tt_nc_g4,   ptr->selection_num_nc, ppt->has_nc_gr4))
          tp_of_tt[index_md][index_tt]=ppt->index_tp_phi_plus_psi;
        if (_index_tt_in_range_(ptr->index_tt_nc_g5,   ptr->selection_num_nc, ppt->has_nc_gr5))
          tp_of_tt[index_md][index_tt]=ppt->index_tp_phi_plus_psi;

        if (_index_tt_in_range_(ptr->index_tt_lensing,  ptr->selection_num_wl,ppt->has_lensing_shear))
          tp_of_tt[index_md][index_tt]=ppt->index_tp_phi_plus_psi;
        if (_index_tt_in_range_(ptr->index_tt_alignment,ptr->selection_num_wl,ppt->has_intrinsic_alignment))
          tp_of_tt[index_md][index_tt]=ppt->index_tp_delta_m;
      }

      if (_vectors_) {

        if ((ppt->has_cl_cmb_temperature == _TRUE_) && (index_tt == ptr->index_tt_t1))
          tp_of_tt[index_md][index_tt]=ppt->index_tp_t1;

        if ((ppt->has_cl_cmb_temperature == _TRUE_) && (index_tt == ptr->index_tt_t2))
          tp_of_tt[index_md][index_tt]=ppt->index_tp_t2;

        if ((ppt->has_cl_cmb_polarization == _TRUE_) && (index_tt == ptr->index_tt_e))
          tp_of_tt[index_md][index_tt]=ppt->index_tp_p;

        if ((ppt->has_cl_cmb_polarization == _TRUE_) && (index_tt == ptr->index_tt_b))
          tp_of_tt[index_md][index_tt]=ppt->index_tp_p;
      }

      if (_tensors_) {

        if ((ppt->has_cl_cmb_temperature == _TRUE_) && (index_tt == ptr->index_tt_t2))
          tp_of_tt[index_md][index_tt]=ppt->index_tp_t2;

        if ((ppt->has_cl_cmb_polarization == _TRUE_) && (index_tt == ptr->index_tt_e))
          tp_of_tt[index_md][index_tt]=ppt->index_tp_p;

        if ((ppt->has_cl_cmb_polarization == _TRUE_) && (index_tt == ptr->index_tt_b))
          tp_of_tt[index_md][index_tt]=ppt->index_tp_p;
      }
    }
  }

  return _SUCCESS_;

}

int transfer_free_source_correspondence(
                                        struct transfers * ptr,
                                        int ** tp_of_tt
                                        ) {

  int index_md;

  for (index_md = 0; index_md < ptr->md_size; index_md++) {
    free(tp_of_tt[index_md]);
  }
  free(tp_of_tt);

  return _SUCCESS_;

}

int transfer_source_tau_size_max(
                                 struct precision * ppr,
                                 struct background * pba,
                                 struct perturbs * ppt,
                                 struct transfers * ptr,
                                 double tau_rec,
                                 double tau0,
                                 int * tau_size_max
                                 ) {

  int index_md;
  int index_tt;
  int tau_size_tt=0;

  *tau_size_max = 0;

  for (index_md = 0; index_md < ptr->md_size; index_md++) {

    for (index_tt = 0; index_tt < ptr->tt_size[index_md]; index_tt++) {

      class_call(transfer_source_tau_size(ppr,
                                          pba,
                                          ppt,
                                          ptr,
                                          tau_rec,
                                          tau0,
                                          index_md,
                                          index_tt,
                                          &tau_size_tt),
                 ptr->error_message,
                 ptr->error_message);

      *tau_size_max = MAX(*tau_size_max,tau_size_tt);
    }
  }

  return _SUCCESS_;
}

/**
 * the code makes a distinction between "perturbation sources"
 * (e.g. gravitational potential) and "transfer sources" (e.g. total
 * density fluctuations, obtained through the Poisson equation, and
 * observed with a given selection function).
 *
 * This routine computes the number of sampled time values for each type
 * of transfer sources.
 *
 * @param ppr                   Input : pointer to precision structure
 * @param pba                   Input : pointer to background structure
 * @param ppt                   Input : pointer to perturbation structure
 * @param ptr                   Input : pointer to transfers structure
 * @param tau_rec               Input : recombination time
 * @param tau0                  Input : time today
 * @param index_md              Input : index of the mode (scalar, tensor)
 * @param index_tt              Input : index of transfer type
 * @param tau_size              Output: pointer to number of smapled times
 * @return the error status
 */

int transfer_source_tau_size(
                             struct precision * ppr,
                             struct background * pba,
                             struct perturbs * ppt,
                             struct transfers * ptr,
                             double tau_rec,
                             double tau0,
                             int index_md,
                             int index_tt,
                             int * tau_size) {

  /* values of conformal time */
  double tau_min,tau_mean,tau_max;

  /* minimum value of index_tt */
  int index_tau_min;

  /* value of l at which limber approximation is switched on */
  int l_limber;

  /* current redhsift bin number */
  int bin=0;

  /* scalar mode */
  if (_scalars_) {

    /* scalar temperature */
    if ((ppt->has_cl_cmb_temperature == _TRUE_) &&
        ((index_tt == ptr->index_tt_t0) || (index_tt == ptr->index_tt_t1) || (index_tt == ptr->index_tt_t2)))
      *tau_size = ppt->tau_size;

    /* scalar polarisation */
    if ((ppt->has_cl_cmb_polarization == _TRUE_) && (index_tt == ptr->index_tt_e))
      *tau_size = ppt->tau_size;

    /* cmb lensing potential */
    if ((ppt->has_cl_cmb_lensing_potential == _TRUE_) && (index_tt == ptr->index_tt_lcmb)) {

      /* find times before recombination, that will be thrown away */
      index_tau_min=0;
      while (ppt->tau_sampling[index_tau_min]<=tau_rec) index_tau_min++;

      /* infer number of time steps after removing early times */
      *tau_size = ppt->tau_size-index_tau_min;
    }

    /* density Cl's */
    if ((_index_tt_in_range_(ptr->index_tt_density, ptr->selection_num_nc, ppt->has_nc_density)) ||
        (_index_tt_in_range_(ptr->index_tt_rsd,     ptr->selection_num_nc, ppt->has_nc_rsd1)) ||
        (_index_tt_in_range_(ptr->index_tt_d0,      ptr->selection_num_nc, ppt->has_nc_rsd2)) ||
        (_index_tt_in_range_(ptr->index_tt_d1,      ptr->selection_num_nc, ppt->has_nc_rsd3)) ||
        (_index_tt_in_range_(ptr->index_tt_nc_g1,   ptr->selection_num_nc, ppt->has_nc_gr1))  ||
        (_index_tt_in_range_(ptr->index_tt_nc_g2,   ptr->selection_num_nc, ppt->has_nc_gr2))  ||
        (_index_tt_in_range_(ptr->index_tt_nc_g3,   ptr->selection_num_nc, ppt->has_nc_gr3))
        ) {

      /* bin number associated to particular redshift bin and selection function */
      if (_index_tt_in_range_(ptr->index_tt_density, ptr->selection_num_nc, ppt->has_nc_density))
        bin = index_tt - ptr->index_tt_density;

      if (_index_tt_in_range_(ptr->index_tt_rsd,     ptr->selection_num_nc, ppt->has_nc_rsd1))
        bin = index_tt - ptr->index_tt_rsd;

      if (_index_tt_in_range_(ptr->index_tt_d0,      ptr->selection_num_nc, ppt->has_nc_rsd2))
        bin = index_tt - ptr->index_tt_d0;

      if (_index_tt_in_range_(ptr->index_tt_d1,      ptr->selection_num_nc, ppt->has_nc_rsd3))
        bin = index_tt - ptr->index_tt_d1;

      if (_index_tt_in_range_(ptr->index_tt_nc_g1,   ptr->selection_num_nc, ppt->has_nc_gr1))
        bin = index_tt - ptr->index_tt_nc_g1;

      if (_index_tt_in_range_(ptr->index_tt_nc_g2,   ptr->selection_num_nc, ppt->has_nc_gr2))
        bin = index_tt - ptr->index_tt_nc_g2;

      if (_index_tt_in_range_(ptr->index_tt_nc_g3,   ptr->selection_num_nc, ppt->has_nc_gr3))
        bin = index_tt - ptr->index_tt_nc_g3;

      /* time interval for this bin */
      class_call(transfer_selection_times_nc(ppr,
					     pba,
					     ppt,
					     ptr,
					     bin,
					     &tau_min,
					     &tau_mean,
					     &tau_max),
                 ptr->error_message,
                 ptr->error_message);

      /* case selection=dirac */
      if (tau_min == tau_max) {
        *tau_size = 1;
      }
      /* other cases (gaussian, top-hat...) */
      else {

        /* check that selection function well sampled */
        *tau_size = (int)ppr->selection_sampling;

        /* value of l at which the code switches to Limber approximation
           (necessary for next step) */
        l_limber=ppr->l_switch_limber_for_cl_density;

        /* check that bessel well sampled, if not define finer sampling
           overwriting the previous one.
           One Bessel oscillations corresponds to [Delta tau]=2pi/k.
           This is minimal for largest relevant k_max,
           namely k_max=l_limber/(tau0-tau_mean).
           We need to cut the interval (tau_max-tau_min) in pieces of size
           [Delta tau]=2pi/k_max. This gives the number below.
        */
        *tau_size=MAX(*tau_size,(int)((tau_max-tau_min)/((tau0-tau_mean)/l_limber))*ppr->selection_sampling_bessel);
      }
    }

    /* galaxy lensing Cl's, differs from density Cl's since the source
       function will spread from the selection function region up to
       tau0 */
    if ((_index_tt_in_range_(ptr->index_tt_nc_lens, ptr->selection_num_nc, ppt->has_nc_lens)) ||
        (_index_tt_in_range_(ptr->index_tt_nc_g4, ptr->selection_num_nc, ppt->has_nc_gr4)) ||
        (_index_tt_in_range_(ptr->index_tt_nc_g5, ptr->selection_num_nc, ppt->has_nc_gr5))
        ) {

      /* bin number associated to particular redshift bin and selection function */
      if (_index_tt_in_range_(ptr->index_tt_nc_lens, ptr->selection_num_nc, ppt->has_nc_lens))
        bin = index_tt - ptr->index_tt_nc_lens;

      if (_index_tt_in_range_(ptr->index_tt_nc_g4,   ptr->selection_num_nc, ppt->has_nc_gr4))
        bin = index_tt - ptr->index_tt_nc_g4;

      if (_index_tt_in_range_(ptr->index_tt_nc_g5,   ptr->selection_num_nc, ppt->has_nc_gr5))
        bin = index_tt - ptr->index_tt_nc_g5;

      /* time interval for this bin */
      class_call(transfer_selection_times_nc(ppr,
					     pba,
					     ppt,
					     ptr,
					     bin,
					     &tau_min,
					     &tau_mean,
					     &tau_max),
                 ptr->error_message,
                 ptr->error_message);

      /* check that selection function well sampled */
      *tau_size = (int)ppr->selection_sampling;

      /* value of l at which the code switches to Limber approximation
         (necessary for next step) */
      l_limber=ppr->l_switch_limber_for_cl_density;

      /* check that bessel well sampled, if not define finer sampling
         overwriting the previous one.
         One Bessel oscillations corresponds to [Delta tau]=2pi/k.
         This is minimal for largest relevant k_max,
         namely k_max=l_limber/((tau0-tau_mean)/2).
         We need to cut the interval (tau_0-tau_min) in pieces of size
         [Delta tau]=2pi/k_max. This gives the number below.
      */
      *tau_size=MAX(*tau_size,(int)((tau0-tau_min)/((tau0-tau_mean)/2./l_limber))*ppr->selection_sampling_bessel);
    }

    /* alignment Cl's */
    if ((_index_tt_in_range_(ptr->index_tt_alignment,ptr->selection_num_wl, ppt->has_intrinsic_alignment))) {
      bin = index_tt - ptr->index_tt_alignment;

      /* time interval for this bin */
      class_call(transfer_selection_times_wl(ppr,
					     pba,
					     ppt,
					     ptr,
					     bin,
					     &tau_min,
					     &tau_mean,
					     &tau_max),
                 ptr->error_message,
                 ptr->error_message);

      /* case selection=dirac */
      if (tau_min == tau_max) {
        *tau_size = 1;
      }
      /* other cases (gaussian, top-hat...) */
      else {

        /* check that selection function well sampled */
        *tau_size = (int)ppr->selection_sampling;

        /* value of l at which the code switches to Limber approximation
           (necessary for next step) */
        l_limber=ppr->l_switch_limber_for_cl_lensing;

        /* check that bessel well sampled, if not define finer sampling
           overwriting the previous one.
           One Bessel oscillations corresponds to [Delta tau]=2pi/k.
           This is minimal for largest relevant k_max,
           namely k_max=l_limber/(tau0-tau_mean).
           We need to cut the interval (tau_max-tau_min) in pieces of size
           [Delta tau]=2pi/k_max. This gives the number below.
        */
        *tau_size=MAX(*tau_size,(int)((tau_max-tau_min)/((tau0-tau_mean)/l_limber))*ppr->selection_sampling_bessel);
      }
    }

    /* lensing shear Cl's, differs from alignment Cl's since the source
       function will spread from the selection function region up to
       tau0 */
    if ((_index_tt_in_range_(ptr->index_tt_lensing,ptr->selection_num_wl, ppt->has_lensing_shear))) {
      bin = index_tt - ptr->index_tt_lensing;

      /* time interval for this bin */
      class_call(transfer_selection_times_wl(ppr,
					     pba,
					     ppt,
					     ptr,
					     bin,
					     &tau_min,
					     &tau_mean,
					     &tau_max),
                 ptr->error_message,
                 ptr->error_message);

      /* check that selection function well sampled */
      *tau_size = (int)ppr->selection_sampling;

      /* value of l at which the code switches to Limber approximation
         (necessary for next step) */
      l_limber=ppr->l_switch_limber_for_cl_lensing;

      /* check that bessel well sampled, if not define finer sampling
         overwriting the previous one.
         One Bessel oscillations corresponds to [Delta tau]=2pi/k.
         This is minimal for largest relevant k_max,
         namely k_max=l_limber/((tau0-tau_mean)/2).
         We need to cut the interval (tau_0-tau_min) in pieces of size
         [Delta tau]=2pi/k_max. This gives the number below.
      */
      *tau_size=MAX(*tau_size,(int)((tau0-tau_min)/((tau0-tau_mean)/2./l_limber))*ppr->selection_sampling_bessel);
    }
  }

  /* tensor mode */
  if (_tensors_) {

    /* for all tensor types */
    *tau_size = ppt->tau_size;
  }

  return _SUCCESS_;
}

int transfer_compute_for_each_q(
                                struct precision * ppr,
                                struct background * pba,
                                struct perturbs * ppt,
                                struct transfers * ptr,
                                int ** tp_of_tt,
                                int index_q,
                                int tau_size_max,
                                double tau_rec,
                                double *** pert_sources,
                                double *** pert_sources_spline,
                                struct transfer_workspace * ptw
                                ) {

  /** Summary: */

  /** - define local variables */

  /* running index for modes */
  int index_md;
  /* running index for initial conditions */
  int index_ic;
  /* running index for transfer types */
  int index_tt;
  /* running index for multipoles */
  int index_l;

  /* we deal with workspaces, i.e. with contiguous memory zones (one
     per thread) containing various fields used by the integration
     routine */

  /* - first workspace field: perturbation source interpolated from perturbation stucture */
  double * interpolated_sources;

  /* - second workspace field: list of tau0-tau values, tau0_minus_tau[index_tau] */
  double * tau0_minus_tau;

  /* - third workspace field: list of trapezoidal weights for integration over tau */
  double * w_trapz;

  /* - fourth workspace field, containing just a double: number of time values */
  int * tau_size;

  /* - fifth workspace field, identical to above interpolated sources:
     sources[index_tau] */
  double * sources;

  /** - for a given l, maximum value of k such that we can convolve
      the source with Bessel functions j_l(x) without reaching x_max */
  double q_max_bessel;

  /* a value of index_type */
  int previous_type;

  double l;

  short neglect;

  radial_function_type radial_type;

  /** store the sources in the workspace and define all
      fields in this workspace */
  interpolated_sources = ptw->interpolated_sources;
  tau0_minus_tau = ptw->tau0_minus_tau;
  w_trapz  = ptw->w_trapz;
  tau_size = &(ptw->tau_size);
  sources = ptw->sources;

  /** - loop over all modes. For each mode: */

  for (index_md = 0; index_md < ptr->md_size; index_md++) {

    /* if we reached q_max for this mode, there is nothing to be done */

    if (ptr->k[index_md][index_q] <= ppt->k[index_md][ppt->k_size_cl[index_md]-1]) {

      /** - loop over initial conditions. For each of them: */

      for (index_ic = 0; index_ic < ppt->ic_size[index_md]; index_ic++) {

        /* initialize the previous type index */
        previous_type=-1;

        /** - loop over types. For each of them: */

        for (index_tt = 0; index_tt < ptr->tt_size[index_md]; index_tt++) {

          /** check if we must now deal with a new source with a
              new index ppt->index_type. If yes, interpolate it at the
              right values of k. */

          if (tp_of_tt[index_md][index_tt] != previous_type) {

            class_call(transfer_interpolate_sources(ppt,
                                                    ptr,
                                                    index_q,
                                                    index_md,
                                                    index_ic,
                                                    tp_of_tt[index_md][index_tt],
                                                    pert_sources[index_md][index_ic * ppt->tp_size[index_md] + tp_of_tt[index_md][index_tt]],
                                                    pert_sources_spline[index_md][index_ic * ppt->tp_size[index_md] + tp_of_tt[index_md][index_tt]],
                                                    interpolated_sources),
                       ptr->error_message,
                       ptr->error_message);
          }

          previous_type = tp_of_tt[index_md][index_tt];

          /* the code makes a distinction between "perturbation
             sources" (e.g. gravitational potential) and "transfer
             sources" (e.g. total density fluctuations, obtained
             through the Poisson equation, and observed with a given
             selection function).

             The next routine computes the transfer source given the
             interpolated perturbation source, and copies it in the
             workspace. */

          class_call(transfer_sources(ppr,
                                      pba,
                                      ppt,
                                      ptr,
                                      interpolated_sources,
                                      tau_rec,
                                      index_q,
                                      index_md,
                                      index_tt,
                                      sources,
                                      tau0_minus_tau,
                                      w_trapz,
                                      tau_size),
                     ptr->error_message,
                     ptr->error_message);

          /* now that the array of times tau0_minus_tau is known, we can
             infer the arry of radial coordinates r(tau0_minus_tau) as well as a
             few other quantities related by trigonometric functions */

          class_call(transfer_radial_coordinates(ptr,ptw,index_md,index_q),
                     ptr->error_message,
                     ptr->error_message);

          /** Select radial function type: */
          class_call(transfer_select_radial_function(
                                                     ppt,
                                                     ptr,
                                                     index_md,
                                                     index_tt,
                                                     &radial_type),
                     ptr->error_message,
                     ptr->error_message);

          for (index_l = 0; index_l < ptr->l_size[index_md]; index_l++) {

            l = (double)ptr->l[index_l];

            /* neglect transfer function when l is much smaller than k*tau0 */
            class_call(transfer_can_be_neglected(ppr,
                                                 ppt,
                                                 ptr,
                                                 index_md,
                                                 index_ic,
                                                 index_tt,
                                                 (pba->conformal_age-tau_rec)*ptr->angular_rescaling,
                                                 ptr->q[index_q],
                                                 l,
                                                 &neglect),
                       ptr->error_message,
                       ptr->error_message);

            /* for K>0 (closed), transfer functions only defined for l<nu */
            if ((ptw->sgnK == 1) && (ptr->l[index_l] >= (int)(ptr->q[index_q]/sqrt(ptw->K)+0.2))) {
              neglect = _TRUE_;
            }
            /* This would maybe go into transfer_can_be_neglected later: */
            if ((ptw->sgnK != 0) && (index_l>=ptw->HIS.l_size) && (index_q < ptr->index_q_flat_approximation)) {
              neglect = _TRUE_;
            }
            if (neglect == _TRUE_) {

              ptr->transfer[index_md][((index_ic * ptr->tt_size[index_md] + index_tt)
                                       * ptr->l_size[index_md] + index_l)
                                      * ptr->q_size + index_q] = 0.;
            }
            else {

              /* for a given l, maximum value of k such that we can
                 convolve the source with Bessel functions j_l(x)
                 without reaching x_max (this is relevant in the flat
                 case when the bessels are compiuted with the old bessel
                 module. otherwise this condition is guaranteed by the
                 choice of proper xmax when computing bessels) */
              if (ptw->sgnK == 0) {
                q_max_bessel = ptw->pBIS->x[ptw->pBIS->x_size-1]/tau0_minus_tau[0];
              }
              else {
                q_max_bessel = ptr->q[ptr->q_size-1];
              }

              /* neglect late time CMB sources when l is above threshold */
              class_call(transfer_late_source_can_be_neglected(ppr,
                                                               ppt,
                                                               ptr,
                                                               index_md,
                                                               index_tt,
                                                               l,
                                                               &(ptw->neglect_late_source)),
                         ptr->error_message,
                         ptr->error_message);

              /* compute the transfer function for this l */
              class_call(transfer_compute_for_each_l(
                                                     ptw,
                                                     ppr,
                                                     ppt,
                                                     ptr,
                                                     index_q,
                                                     index_md,
                                                     index_ic,
                                                     index_tt,
                                                     index_l,
                                                     l,
                                                     q_max_bessel,
                                                     radial_type
                                                     ),
                         ptr->error_message,
                         ptr->error_message);
            }

          } /* end of loop over l */

        } /* end of loop over type */

      } /* end of loop over initial condition */

    }

    else {

      for (index_ic = 0; index_ic < ppt->ic_size[index_md]; index_ic++) {
        for (index_tt = 0; index_tt < ptr->tt_size[index_md]; index_tt++) {
          for (index_l = 0; index_l < ptr->l_size[index_md]; index_l++) {

            ptr->transfer[index_md][((index_ic * ptr->tt_size[index_md] + index_tt)
                                     * ptr->l_size[index_md] + index_l)
                                    * ptr->q_size + index_q] = 0.;
          }
        }
      }
    }

  } /* end of loop over mode */

  return _SUCCESS_;

}

int transfer_gather_for_each_q(
                                struct precision * ppr,
                                struct background * pba,
                                struct perturbs * ppt,
                                struct transfers * ptr,
                                int index_q)
{

  /** Summary: */

  /** - define local variables */

  /* running index for modes */
  int index_md;
  /* running index for initial conditions */
  int index_ic;
  /* running index for transfer types */
  int index_tt;
  /* running index for multipoles */
  int index_l;
  /* running index for redshift bins */
  int index_d1;
  /* ell for each index_l */
  double l;

  /** - loop over all modes. For each mode: */
  for (index_md = 0; index_md < ptr->md_size; index_md++) {
    if (ptr->k[index_md][index_q] <= ppt->k[index_md][ppt->k_size_cl[index_md]-1]) {
      for (index_ic = 0; index_ic < ppt->ic_size[index_md]; index_ic++) {
	for (index_l = 0; index_l < ptr->l_size[index_md]; index_l++) {
	  l = (double)ptr->l[index_l];
	  
	  if(_scalars_) {
	    if(ppt->has_cl_cmb_temperature==_TRUE_) { //T
	      ptr->transfer_all[index_md][((index_ic*ptr->tt_size_all[index_md]+ptr->index_tt_tall)
					   *ptr->l_size[index_md]+index_l)*ptr->q_size+index_q]=
		ptr->transfer[index_md][((index_ic*ptr->tt_size[index_md]+ptr->index_tt_t0)
					 *ptr->l_size[index_md]+index_l)*ptr->q_size+index_q]+
		ptr->transfer[index_md][((index_ic*ptr->tt_size[index_md]+ptr->index_tt_t1)
					 *ptr->l_size[index_md]+index_l)*ptr->q_size+index_q]+
		ptr->transfer[index_md][((index_ic*ptr->tt_size[index_md]+ptr->index_tt_t2)
					 *ptr->l_size[index_md]+index_l)*ptr->q_size+index_q];
	    }
	    if(ppt->has_cl_cmb_polarization==_TRUE_) { //E
	      ptr->transfer_all[index_md][((index_ic*ptr->tt_size_all[index_md]+ptr->index_tt_eall)
					   *ptr->l_size[index_md]+index_l)*ptr->q_size+index_q]=
		ptr->transfer[index_md][((index_ic*ptr->tt_size[index_md]+ptr->index_tt_e)
					 *ptr->l_size[index_md]+index_l)*ptr->q_size+index_q];
	    }
	    if(ppt->has_cl_cmb_lensing_potential==_TRUE_) { //P
	      ptr->transfer_all[index_md][((index_ic*ptr->tt_size_all[index_md]+ptr->index_tt_pall)
					   *ptr->l_size[index_md]+index_l)*ptr->q_size+index_q]=
		ptr->transfer[index_md][((index_ic*ptr->tt_size[index_md]+ptr->index_tt_lcmb)
					 *ptr->l_size[index_md]+index_l)*ptr->q_size+index_q];
	    }
	    if(ppt->has_cl_number_count==_TRUE_) { //NC
	      for(index_d1=0;index_d1<ptr->selection_num_nc;index_d1++) {
		double transf=0;
		if(ppt->has_nc_density==_TRUE_) {
		  transf+=
		    ptr->transfer[index_md][((index_ic*ptr->tt_size[index_md]+ptr->index_tt_density+index_d1)
					     *ptr->l_size[index_md]+index_l)*ptr->q_size+index_q];
		}
		if(ppt->has_nc_rsd1==_TRUE_) {
		  transf+=
		    ptr->transfer[index_md][((index_ic*ptr->tt_size[index_md]+ptr->index_tt_rsd+index_d1)
					     *ptr->l_size[index_md]+index_l)*ptr->q_size+index_q];
		}
		if(ppt->has_nc_rsd2==_TRUE_) {
		  transf+=
		    ptr->transfer[index_md][((index_ic*ptr->tt_size[index_md]+ptr->index_tt_d0+index_d1)
					     *ptr->l_size[index_md]+index_l)*ptr->q_size+index_q];
		}
		if(ppt->has_nc_rsd3==_TRUE_) {
		  transf+=
		    ptr->transfer[index_md][((index_ic*ptr->tt_size[index_md]+ptr->index_tt_d1+index_d1)
					     *ptr->l_size[index_md]+index_l)*ptr->q_size+index_q];
		}
		if(ppt->has_nc_lens==_TRUE_) {
		  transf+=
		    l*(l+1)*ptr->transfer[index_md][((index_ic*ptr->tt_size[index_md]+ptr->index_tt_nc_lens+index_d1)
						     *ptr->l_size[index_md]+index_l)*ptr->q_size+index_q];
		}
		if(ppt->has_nc_gr1==_TRUE_) {
		  transf+=
		    ptr->transfer[index_md][((index_ic*ptr->tt_size[index_md]+ptr->index_tt_nc_g1+index_d1)
						     *ptr->l_size[index_md]+index_l)*ptr->q_size+index_q];
		}
		if(ppt->has_nc_gr2==_TRUE_) {
		  transf+=
		    ptr->transfer[index_md][((index_ic*ptr->tt_size[index_md]+ptr->index_tt_nc_g2+index_d1)
						     *ptr->l_size[index_md]+index_l)*ptr->q_size+index_q];
		}
		if(ppt->has_nc_gr3==_TRUE_) {
		  transf+=
		    ptr->transfer[index_md][((index_ic*ptr->tt_size[index_md]+ptr->index_tt_nc_g3+index_d1)
						     *ptr->l_size[index_md]+index_l)*ptr->q_size+index_q];
		}
		if(ppt->has_nc_gr4==_TRUE_) {
		  transf+=
		    ptr->transfer[index_md][((index_ic*ptr->tt_size[index_md]+ptr->index_tt_nc_g4+index_d1)
						     *ptr->l_size[index_md]+index_l)*ptr->q_size+index_q];
		}
		if(ppt->has_nc_gr5==_TRUE_) {
		  transf+=
		    ptr->transfer[index_md][((index_ic*ptr->tt_size[index_md]+ptr->index_tt_nc_g5+index_d1)
						     *ptr->l_size[index_md]+index_l)*ptr->q_size+index_q];
		}
		ptr->transfer_all[index_md][((index_ic*ptr->tt_size_all[index_md]+ptr->index_tt_ncall+index_d1)
					     *ptr->l_size[index_md]+index_l)*ptr->q_size+index_q]=transf;
	      }
	    }
	    if(ppt->has_weak_lensing==_TRUE_) { //WL
	      for(index_d1=0;index_d1<ptr->selection_num_wl;index_d1++) {
		double transf=0;
		if(ppt->has_intrinsic_alignment==_TRUE_) {
		  transf+=
		    ptr->transfer[index_md][((index_ic*ptr->tt_size[index_md]+ptr->index_tt_alignment+index_d1)
					     *ptr->l_size[index_md]+index_l)*ptr->q_size+index_q];
		}
		if(ppt->has_lensing_shear==_TRUE_) {
		  transf+=
		    ptr->transfer[index_md][((index_ic*ptr->tt_size[index_md]+ptr->index_tt_lensing+index_d1)
					     *ptr->l_size[index_md]+index_l)*ptr->q_size+index_q];
		}
		transf*=sqrt((l+2.)*(l+1.)*l*(l-1.));
		ptr->transfer_all[index_md][((index_ic*ptr->tt_size_all[index_md]+ptr->index_tt_wlall+index_d1)
					     *ptr->l_size[index_md]+index_l)*ptr->q_size+index_q]=transf;
	      }
	    }
	  }
	  if(_vectors_) {
	    if(ppt->has_cl_cmb_temperature==_TRUE_) { //T
	      ptr->transfer_all[index_md][((index_ic*ptr->tt_size_all[index_md]+ptr->index_tt_tall)
					   *ptr->l_size[index_md]+index_l)*ptr->q_size+index_q]=
		ptr->transfer[index_md][((index_ic*ptr->tt_size[index_md]+ptr->index_tt_t1)
					 *ptr->l_size[index_md]+index_l)*ptr->q_size+index_q]+
		ptr->transfer[index_md][((index_ic*ptr->tt_size[index_md]+ptr->index_tt_t2)
					 *ptr->l_size[index_md]+index_l)*ptr->q_size+index_q];
	    }
	    if(ppt->has_cl_cmb_polarization==_TRUE_) { //E,B
	      ptr->transfer_all[index_md][((index_ic*ptr->tt_size_all[index_md]+ptr->index_tt_eall)
					   *ptr->l_size[index_md]+index_l)*ptr->q_size+index_q]=
		ptr->transfer[index_md][((index_ic*ptr->tt_size[index_md]+ptr->index_tt_e)
					 *ptr->l_size[index_md]+index_l)*ptr->q_size+index_q];
	      ptr->transfer_all[index_md][((index_ic*ptr->tt_size_all[index_md]+ptr->index_tt_ball)
					   *ptr->l_size[index_md]+index_l)*ptr->q_size+index_q]=
		ptr->transfer[index_md][((index_ic*ptr->tt_size[index_md]+ptr->index_tt_b)
					 *ptr->l_size[index_md]+index_l)*ptr->q_size+index_q];
	    }
	  }
	  if(_tensors_) {
	    if(ppt->has_cl_cmb_temperature==_TRUE_) { //T
	      ptr->transfer_all[index_md][((index_ic*ptr->tt_size_all[index_md]+ptr->index_tt_tall)
					   *ptr->l_size[index_md]+index_l)*ptr->q_size+index_q]=
		ptr->transfer[index_md][((index_ic*ptr->tt_size[index_md]+ptr->index_tt_t2)
					 *ptr->l_size[index_md]+index_l)*ptr->q_size+index_q];
	    }
	    if(ppt->has_cl_cmb_polarization==_TRUE_) { //E,B
	      ptr->transfer_all[index_md][((index_ic*ptr->tt_size_all[index_md]+ptr->index_tt_eall)
					   *ptr->l_size[index_md]+index_l)*ptr->q_size+index_q]=
		ptr->transfer[index_md][((index_ic*ptr->tt_size[index_md]+ptr->index_tt_e)
					 *ptr->l_size[index_md]+index_l)*ptr->q_size+index_q];
	      ptr->transfer_all[index_md][((index_ic*ptr->tt_size_all[index_md]+ptr->index_tt_ball)
					   *ptr->l_size[index_md]+index_l)*ptr->q_size+index_q]=
		ptr->transfer[index_md][((index_ic*ptr->tt_size[index_md]+ptr->index_tt_b)
					 *ptr->l_size[index_md]+index_l)*ptr->q_size+index_q];
	    }
	  } 

	} /* end of loop over l */
      } /* end of loop over initial condition */
    }
    else {
      for (index_ic = 0; index_ic < ppt->ic_size[index_md]; index_ic++) {
        for (index_tt = 0; index_tt < ptr->tt_size_all[index_md]; index_tt++) {
          for (index_l = 0; index_l < ptr->l_size[index_md]; index_l++) {
            ptr->transfer_all[index_md][((index_ic*ptr->tt_size_all[index_md]+index_tt)
					 *ptr->l_size[index_md]+index_l)
					*ptr->q_size+index_q]=0.;
          }
        }
      }
    }
  } /* end of loop over mode */

  return _SUCCESS_;

}

int transfer_radial_coordinates(
                                struct transfers * ptr,
                                struct transfer_workspace * ptw,
                                int index_md,
                                int index_q
                                ) {

  int index_tau;
  double sqrt_absK=0.;

  switch (ptw->sgnK){
  case 1:
    sqrt_absK = sqrt(ptw->K);
    for (index_tau=0; index_tau < ptw->tau_size; index_tau++) {
      ptw->chi[index_tau] = sqrt_absK*ptw->tau0_minus_tau[index_tau];
      ptw->cscKgen[index_tau] = sqrt_absK/ptr->k[index_md][index_q]/sin(ptw->chi[index_tau]);
      ptw->cotKgen[index_tau] = ptw->cscKgen[index_tau]*cos(ptw->chi[index_tau]);
    }
    break;
  case 0:
    for (index_tau=0; index_tau < ptw->tau_size; index_tau++) {
      ptw->chi[index_tau] = ptr->k[index_md][index_q] * ptw->tau0_minus_tau[index_tau];
      ptw->cscKgen[index_tau] = 1.0/ptw->chi[index_tau];
      ptw->cotKgen[index_tau] = 1.0/ptw->chi[index_tau];
    }
    break;
  case -1:
    sqrt_absK = sqrt(-ptw->K);
    for (index_tau=0; index_tau < ptw->tau_size; index_tau++) {
      ptw->chi[index_tau] = sqrt_absK*ptw->tau0_minus_tau[index_tau];
      ptw->cscKgen[index_tau] = sqrt_absK/ptr->k[index_md][index_q]/sinh(ptw->chi[index_tau]);
      ptw->cotKgen[index_tau] = ptw->cscKgen[index_tau]*cosh(ptw->chi[index_tau]);
    }
    break;
  }

  return _SUCCESS_;
}


/**
 * This routine interpolates sources \f$ S(k, \tau) \f$ for each mode,
 * initial condition and type (of perturbation module), to get them at
 * the right values of k, using the spline interpolation method.
 *
 * @param ppt                   Input : pointer to perturbation structure
 * @param ptr                   Input : pointer to transfers structure
 * @param index_md              Input : index of mode
 * @param index_ic              Input : index of initial condition
 * @param index_type            Input : index of type of source (in perturbation module)
 * @param pert_source           Input : array of sources
 * @param pert_source_spline    Input : array of second derivative of sources
 * @param interpolated_sources  Output: array of interpolated sources (filled here but allocated in transfer_init() to avoid numerous reallocation)
 * @return the error status
 */

int transfer_interpolate_sources(
                                 struct perturbs * ppt,
                                 struct transfers * ptr,
                                 int index_q,
                                 int index_md,
                                 int index_ic,
                                 int index_type,
                                 double * pert_source,       /* array with argument pert_source[index_tau*ppt->k_size[index_md]+index_k] (must be allocated) */
                                 double * pert_source_spline, /* array with argument pert_source_spline[index_tau*ppt->k_size[index_md]+index_k] (must be allocated) */
                                 double * interpolated_sources /* array with argument interpolated_sources[index_q*ppt->tau_size+index_tau] (must be allocated) */
                                 ) {

  /** Summary: */

  /** - define local variables */

  /* index running on k values in the original source array */
  int index_k;

  /* index running on time */
  int index_tau;

  /* variables used for spline interpolation algorithm */
  double h, a, b;

  /** - interpolate at each k value using the usual
      spline interpolation algorithm. */

  index_k = 0;
  h = ppt->k[index_md][index_k+1] - ppt->k[index_md][index_k];

  while (((index_k+1) < ppt->k_size[index_md]) &&
         (ppt->k[index_md][index_k+1] <
          ptr->k[index_md][index_q])) {
    index_k++;
    h = ppt->k[index_md][index_k+1] - ppt->k[index_md][index_k];
  }

  class_test(h==0.,
             ptr->error_message,
             "stop to avoid division by zero");

  b = (ptr->k[index_md][index_q] - ppt->k[index_md][index_k])/h;
  a = 1.-b;

  for (index_tau = 0; index_tau < ppt->tau_size; index_tau++) {

    interpolated_sources[index_tau] =
      a * pert_source[index_tau*ppt->k_size[index_md]+index_k]
      + b * pert_source[index_tau*ppt->k_size[index_md]+index_k+1]
      + ((a*a*a-a) * pert_source_spline[index_tau*ppt->k_size[index_md]+index_k]
         +(b*b*b-b) * pert_source_spline[index_tau*ppt->k_size[index_md]+index_k+1])*h*h/6.0;

  }

  return _SUCCESS_;

}

/**
 * the code makes a distinction between "perturbation sources"
 * (e.g. gravitational potential) and "transfer sources" (e.g. total
 * density fluctuations, obtained through the Poisson equation, and
 * observed with a given selection function).
 *
 * This routine computes the transfer source given the interpolated
 * perturbation source, and copies it in the workspace.
 *
 * @param ppr                   Input : pointer to precision structure
 * @param pba                   Input : pointer to background structure
 * @param ppt                   Input : pointer to perturbation structure
 * @param ptr                   Input : pointer to transfers structure
 * @param interpolated_sources  Input : interpolated perturbation source
 * @param tau_rec               Input : recombination time
 * @param index_md            Input : index of mode
 * @param index_tt              Input : index of type of (transfer) source
 * @param sources               Output: transfer source
 * @param tau0_minus_tau        Output: values of (tau0-tau) at which source are sample
 * @param w_trapz               Output: trapezoidal weights for integration over tau
 * @param tau_size_double       Output: pointer to size of previous two arrays, converted to double
 * @return the error status
 */

int transfer_sources(
                     struct precision * ppr,
                     struct background * pba,
                     struct perturbs * ppt,
                     struct transfers * ptr,
                     double * interpolated_sources,
                     double tau_rec,
                     int index_q,
                     int index_md,
                     int index_tt,
                     double * sources,
                     double * tau0_minus_tau,
                     double * w_trapz,
                     int * tau_size_out
                     )  {

  /** Summary: */

  /** - define local variables */

  /* index running on time */
  int index_tau;

  /* bin for computation of cl_density */
  int bin=0;
  short itr;

  /* number of tau values */
  int tau_size;

  /* minimum tau index kept in transfer sources */
  int index_tau_min;

  /* for calling background_at_eta */
  int last_index;
  double * pvecback = NULL;

  /* conformal time */
  double tau, tau0;

  /* rescaling factor depending on the background at a given time */
  double rescaling=0.;

  /* flag: is there any difference between the perturbation and transfer source? */
  short redefine_source;

  /* array of selection function values at different times */
  double * selection;

  /* clustering, magnification and evolution biases */
  double b_bias=1.,e_bias=0.,s_bias=0.4;

  /* when the selection function is multiplied by a function dNdz */
  double z;

  /* in which cases are perturbation and transfer sources are different?
     I.e., in which case do we need to mutiply the sources by some
     background and/or window function, and eventually to resample it,
     or redfine its time limits? */

  redefine_source = _FALSE_;

  if (_scalars_) {

    /* cmb lensing potential */
    if ((ppt->has_cl_cmb_lensing_potential == _TRUE_) && (index_tt == ptr->index_tt_lcmb))
      redefine_source = _TRUE_;

    /* number count Cl's */
    if ((_index_tt_in_range_(ptr->index_tt_density, ptr->selection_num_nc, ppt->has_nc_density)) ||
        (_index_tt_in_range_(ptr->index_tt_rsd,     ptr->selection_num_nc, ppt->has_nc_rsd1)) ||
        (_index_tt_in_range_(ptr->index_tt_d0,      ptr->selection_num_nc, ppt->has_nc_rsd2)) ||
        (_index_tt_in_range_(ptr->index_tt_d1,      ptr->selection_num_nc, ppt->has_nc_rsd3)) ||
        (_index_tt_in_range_(ptr->index_tt_nc_lens, ptr->selection_num_nc, ppt->has_nc_lens))||
        (_index_tt_in_range_(ptr->index_tt_nc_g1,   ptr->selection_num_nc, ppt->has_nc_gr1))  ||
        (_index_tt_in_range_(ptr->index_tt_nc_g2,   ptr->selection_num_nc, ppt->has_nc_gr2))  ||
        (_index_tt_in_range_(ptr->index_tt_nc_g3,   ptr->selection_num_nc, ppt->has_nc_gr3))  ||
        (_index_tt_in_range_(ptr->index_tt_nc_g4,   ptr->selection_num_nc, ppt->has_nc_gr4))  ||
        (_index_tt_in_range_(ptr->index_tt_nc_g5,   ptr->selection_num_nc, ppt->has_nc_gr5))
        )
      redefine_source = _TRUE_;

    /* galaxy lensing potential */
    if ((_index_tt_in_range_(ptr->index_tt_alignment,ptr->selection_num_wl,ppt->has_intrinsic_alignment)) ||
	(_index_tt_in_range_(ptr->index_tt_lensing  ,ptr->selection_num_wl,ppt->has_lensing_shear))
	)
      redefine_source = _TRUE_;
  }

  /* conformal time today */
  tau0 = pba->conformal_age;

  /* case where we need to redefine by a window function (or any
     function of the background and of k) */
  if (redefine_source == _TRUE_) {

    class_call(transfer_source_tau_size(ppr,
                                        pba,
                                        ppt,
                                        ptr,
                                        tau_rec,
                                        tau0,
                                        index_md,
                                        index_tt,
                                        &tau_size),
               ptr->error_message,
               ptr->error_message);

    if (_scalars_) {

      /* lensing source: throw away times before recombuination, and multiply psi by window function */

      if ((ppt->has_cl_cmb_lensing_potential == _TRUE_) && (index_tt == ptr->index_tt_lcmb)) {

        /* first time step after removing early times */
        index_tau_min =  ppt->tau_size - tau_size;

        /* loop over time and rescale */
        for (index_tau = index_tau_min; index_tau < ppt->tau_size; index_tau++) {

          /* conformal time */
          tau = ppt->tau_sampling[index_tau];

          /* lensing source =  - W(tau) (phi(k,tau) + psi(k,tau)) Heaviside(tau-tau_rec)
             with
             psi,phi = metric perturbation in newtonian gauge (phi+psi = Phi_A-Phi_H of Bardeen)
             W = (tau-tau_rec)/(tau_0-tau)/(tau_0-tau_rec)
             H(x) = Heaviside
             (in tau = tau_0, set source = 0 to avoid division by zero;
             regulated anyway by Bessel).
          */

          if (index_tau == ppt->tau_size-1) {
            rescaling=0.;
          }
          else {
            switch (pba->sgnK){
            case 1:
              rescaling = sqrt(pba->K)
                *sin((tau_rec-tau)*sqrt(pba->K))
                /sin((tau0-tau)*sqrt(pba->K))
                /sin((tau0-tau_rec)*sqrt(pba->K));
              break;
            case 0:
              rescaling = (tau_rec-tau)/(tau0-tau)/(tau0-tau_rec);
              break;
            case -1:
              rescaling = sqrt(-pba->K)
                *sinh((tau_rec-tau)*sqrt(-pba->K))
                /sinh((tau0-tau)*sqrt(-pba->K))
                /sinh((tau0-tau_rec)*sqrt(-pba->K));
              break;
            }
            // Note: until 2.4.3 there was a bug here: the curvature effects had been ommitted.
          }

          /* copy from input array to output array */
          sources[index_tau-index_tau_min] =
            interpolated_sources[index_tau]
            * rescaling
            * ptr->lcmb_rescale
            * pow(ptr->k[index_md][index_q]/ptr->lcmb_pivot,ptr->lcmb_tilt);

          /* store value of (tau0-tau) */
          tau0_minus_tau[index_tau-index_tau_min] = tau0 - tau;

        }

        /* Compute trapezoidal weights for integration over tau */
        class_call(array_trapezoidal_mweights(tau0_minus_tau,
                                              tau_size,
                                              w_trapz,
                                              ptr->error_message),
                   ptr->error_message,
                   ptr->error_message);
      }

      /* density source: redefine the time sampling, multiply by
         coefficient of Poisson equation, and multiply by selection
         function */
      if ((_index_tt_in_range_(ptr->index_tt_density, ptr->selection_num_nc, ppt->has_nc_density)) ||
          (_index_tt_in_range_(ptr->index_tt_rsd,     ptr->selection_num_nc, ppt->has_nc_rsd1)) ||
          (_index_tt_in_range_(ptr->index_tt_d0,      ptr->selection_num_nc, ppt->has_nc_rsd2)) ||
          (_index_tt_in_range_(ptr->index_tt_d1,      ptr->selection_num_nc, ppt->has_nc_rsd3)) ||
          (_index_tt_in_range_(ptr->index_tt_nc_g1,   ptr->selection_num_nc, ppt->has_nc_gr1))  ||
          (_index_tt_in_range_(ptr->index_tt_nc_g2,   ptr->selection_num_nc, ppt->has_nc_gr2))  ||
          (_index_tt_in_range_(ptr->index_tt_nc_g3,   ptr->selection_num_nc, ppt->has_nc_gr3))
          ) {

        /* bin number associated to particular redshift bin and selection function */
        if (_index_tt_in_range_(ptr->index_tt_density, ptr->selection_num_nc, ppt->has_nc_density))
          bin = index_tt - ptr->index_tt_density;

        if (_index_tt_in_range_(ptr->index_tt_rsd,     ptr->selection_num_nc, ppt->has_nc_rsd1))
          bin = index_tt - ptr->index_tt_rsd;

        if (_index_tt_in_range_(ptr->index_tt_d0,      ptr->selection_num_nc, ppt->has_nc_rsd2))
          bin = index_tt - ptr->index_tt_d0;

        if (_index_tt_in_range_(ptr->index_tt_d1,      ptr->selection_num_nc, ppt->has_nc_rsd3))
          bin = index_tt - ptr->index_tt_d1;

        if (_index_tt_in_range_(ptr->index_tt_nc_g1,   ptr->selection_num_nc, ppt->has_nc_gr1))
          bin = index_tt - ptr->index_tt_nc_g1;

        if (_index_tt_in_range_(ptr->index_tt_nc_g2,   ptr->selection_num_nc, ppt->has_nc_gr2))
          bin = index_tt - ptr->index_tt_nc_g2;

        if (_index_tt_in_range_(ptr->index_tt_nc_g3,   ptr->selection_num_nc, ppt->has_nc_gr3))
          bin = index_tt - ptr->index_tt_nc_g3;

        /* allocate temporary arrays for storing sources and for calling background */
        class_alloc(pvecback,pba->bg_size*sizeof(double),ptr->error_message);

        /* redefine the time sampling */
	memcpy(tau0_minus_tau,ptr->glb_selection_0_tau0_minus_tau[bin],ptr->glb_selection_0_size[bin]*sizeof(double));

        /* resample the source at those times */
        class_call(transfer_source_resample(ppr,
                                            pba,
                                            ppt,
                                            ptr,
                                            bin,
                                            tau0_minus_tau,
                                            tau_size,
                                            index_md,
                                            tau0,
                                            interpolated_sources,
                                            sources),
                   ptr->error_message,
                   ptr->error_message);

        /* Compute trapezoidal weights for integration over tau */
        class_call(array_trapezoidal_mweights(tau0_minus_tau,
                                              tau_size,
                                              w_trapz,
                                              ptr->error_message),
                   ptr->error_message,
                   ptr->error_message);

        /* compute values of selection function at sampled values of tau */
	selection=ptr->glb_selection_0[bin];
	
	//Compute tracer type
	itr=ptr->selection_tracer_nc[bin];

        /* loop over time and rescale */
        for (index_tau = 0; index_tau < tau_size; index_tau++) {

          /* conformal time */
	  tau = tau0 - tau0_minus_tau[index_tau];

          /* corresponding background quantities */
          class_call(background_at_tau(pba,
				       tau,
				       pba->long_info,
                                       pba->inter_normal,
                                       &last_index,
                                       pvecback),
                     pba->error_message,
                     ptr->error_message);
	  
	  /* Clustering bias for density term */
          if (_index_tt_in_range_(ptr->index_tt_density, ptr->selection_num_nc, ppt->has_nc_density)) {
	    z = pba->a_today/pvecback[pba->index_bg_a]-1.;
	    class_call(transfer_interpolate_bias_function(ppr,pba,ppt,ptr,
							  z,itr,
							  ptr->bz_z,ptr->bz_bz,ptr->bz_ddbz,ptr->bz_size,
							  &last_index,&b_bias),
		       ptr->error_message,ptr->error_message);
	  }

          /* Magnification bias, used by number counf rsd and number count gravity terms */
          if ((_index_tt_in_range_(ptr->index_tt_d1,      ptr->selection_num_nc, ppt->has_nc_rsd3)) ||
              (_index_tt_in_range_(ptr->index_tt_nc_g2,   ptr->selection_num_nc, ppt->has_nc_gr2))) {
	    z = pba->a_today/pvecback[pba->index_bg_a]-1.;
	    class_call(transfer_interpolate_bias_function(ppr,pba,ppt,ptr,
							  z,itr,
							  ptr->sz_z,ptr->sz_sz,ptr->sz_ddsz,ptr->sz_size,
							  &last_index,&s_bias),
		       ptr->error_message,ptr->error_message);
	  }

          /* Source evolution, used by number counf rsd and number count gravity terms */
          if ((_index_tt_in_range_(ptr->index_tt_d0,      ptr->selection_num_nc, ppt->has_nc_rsd2)) ||
	      (_index_tt_in_range_(ptr->index_tt_d1,      ptr->selection_num_nc, ppt->has_nc_rsd3)) ||
              (_index_tt_in_range_(ptr->index_tt_nc_g2,   ptr->selection_num_nc, ppt->has_nc_gr2))) {
	    z = pba->a_today/pvecback[pba->index_bg_a]-1.;
	    class_call(transfer_interpolate_bias_function(ppr,pba,ppt,ptr,
							  z,itr,
							  ptr->ez_z,ptr->ez_ez,ptr->ez_ddez,ptr->ez_size,
							  &last_index,&e_bias),
		       ptr->error_message,ptr->error_message);
	  }

          /* matter density source =  [- (dz/dtau) W(z)] * delta_m(k,tau)
             = W(tau) delta_m(k,tau)
             with
             delta_m = total matter perturbation (defined in gauge-independent way, see arXiv 1307.1459)
             W(z) = redshift space selection function = dN/dz
             W(tau) = same wrt conformal time = dN/dtau
             (in tau = tau_0, set source = 0 to avoid division by zero;
             regulated anyway by Bessel).
          */

          if (_index_tt_in_range_(ptr->index_tt_density, ptr->selection_num_nc, ppt->has_nc_density)) {
	    double bk;
	    if(ppt->do_f_nl==_TRUE_) {
	      double Dz=pvecback[pba->index_bg_D]*ppt->inv_growth_0;
	      double fnl_kdep=get_f_nl_at_k(ptr->k[index_md][index_q],ppt);
	      bk=b_bias+(b_bias-1)*fnl_kdep/Dz;
	    }
	    else
	      bk=b_bias;
            rescaling = bk*selection[index_tau];
	  }

          /* redhsift space distorsion source = - [- (dz/dtau) W(z)] * (k/H) * theta(k,tau) */
          if (_index_tt_in_range_(ptr->index_tt_rsd,     ptr->selection_num_nc, ppt->has_nc_rsd1))
            rescaling = selection[index_tau]/pvecback[pba->index_bg_H]/pvecback[pba->index_bg_a];

          if (_index_tt_in_range_(ptr->index_tt_d0,      ptr->selection_num_nc, ppt->has_nc_rsd2))
            rescaling = (e_bias-3.)*selection[index_tau]*pvecback[pba->index_bg_H]*pvecback[pba->index_bg_a]
              /ptr->k[index_md][index_q]/ptr->k[index_md][index_q];

          if (_index_tt_in_range_(ptr->index_tt_d1,      ptr->selection_num_nc, ppt->has_nc_rsd3))
            rescaling = selection[index_tau]*(1.
                                              +pvecback[pba->index_bg_H_prime]
                                              /pvecback[pba->index_bg_a]
                                              /pvecback[pba->index_bg_H]
                                              /pvecback[pba->index_bg_H]
                                              +(2.-5.*s_bias)
                                              /tau0_minus_tau[index_tau]
                                              /pvecback[pba->index_bg_a]
                                              /pvecback[pba->index_bg_H]
                                              +5.*s_bias
                                              -e_bias
                                              )/ptr->k[index_md][index_q];

          if (_index_tt_in_range_(ptr->index_tt_nc_g1,   ptr->selection_num_nc, ppt->has_nc_gr1))
            rescaling = selection[index_tau]; //DAM: should there be a 5s-2 factor here??

          if (_index_tt_in_range_(ptr->index_tt_nc_g2,   ptr->selection_num_nc, ppt->has_nc_gr2))
            rescaling = -selection[index_tau]*(3. //DAM: should this be a 2+5s ??, and why the negative sign??
                                               +pvecback[pba->index_bg_H_prime]
                                               /pvecback[pba->index_bg_a]
                                               /pvecback[pba->index_bg_H]
                                               /pvecback[pba->index_bg_H]
                                               +(2.-5.*s_bias)
                                               /tau0_minus_tau[index_tau]
                                               /pvecback[pba->index_bg_a]
                                               /pvecback[pba->index_bg_H]
                                               -e_bias
                                               );

          if (_index_tt_in_range_(ptr->index_tt_nc_g3,   ptr->selection_num_nc, ppt->has_nc_gr3))
            rescaling = selection[index_tau]/pvecback[pba->index_bg_a]/pvecback[pba->index_bg_H];

          sources[index_tau] *= rescaling;
        }

        /* deallocate temporary arrays */
        free(pvecback);
      }

      /* lensing potential: eliminate early times, and multiply by selection
         function */
      if ((_index_tt_in_range_(ptr->index_tt_nc_lens, ptr->selection_num_nc, ppt->has_nc_lens)) ||
          (_index_tt_in_range_(ptr->index_tt_nc_g4,   ptr->selection_num_nc, ppt->has_nc_gr4)) ||
          (_index_tt_in_range_(ptr->index_tt_nc_g5,   ptr->selection_num_nc, ppt->has_nc_gr5))
          ) {

        if (_index_tt_in_range_(ptr->index_tt_nc_lens, ptr->selection_num_nc, ppt->has_nc_lens))
          bin = index_tt - ptr->index_tt_nc_lens;

        if (_index_tt_in_range_(ptr->index_tt_nc_g4,   ptr->selection_num_nc, ppt->has_nc_gr4))
          bin = index_tt - ptr->index_tt_nc_g4;

        if (_index_tt_in_range_(ptr->index_tt_nc_g5,   ptr->selection_num_nc, ppt->has_nc_gr5))
          bin = index_tt - ptr->index_tt_nc_g5;

        /* redefine the time sampling */
	memcpy(tau0_minus_tau,ptr->glb_selection_lns_tau0_minus_tau_nc[bin],
	       ptr->glb_selection_lns_size_nc[bin]*sizeof(double));

        /* resample the source at those times */
        class_call(transfer_source_resample(ppr,
                                            pba,
                                            ppt,
                                            ptr,
                                            bin,
                                            tau0_minus_tau,
                                            tau_size,
                                            index_md,
                                            tau0,
                                            interpolated_sources,
                                            sources),
                   ptr->error_message,
                   ptr->error_message);

        /* Compute trapezoidal weights for integration over tau */
        class_call(array_trapezoidal_mweights(tau0_minus_tau,
                                              tau_size,
                                              w_trapz,
                                              ptr->error_message),
                   ptr->error_message,
                   ptr->error_message);

        /* loop over time and rescale */
        for (index_tau = 0; index_tau < tau_size; index_tau++) {
          /* lensing source =  - W(tau) (phi(k,tau) + psi(k,tau)) Heaviside(tau-tau_rec)
             with
             psi,phi = metric perturbation in newtonian gauge (phi+psi = Phi_A-Phi_H of Bardeen)
             W = (tau-tau_rec)/(tau_0-tau)/(tau_0-tau_rec)
             H(x) = Heaviside
             (in tau = tau_0, set source = 0 to avoid division by zero;
             regulated anyway by Bessel).
          */
	  
	  if (_index_tt_in_range_(ptr->index_tt_nc_lens, ptr->selection_num_nc, ppt->has_nc_lens)) {
	    rescaling=-ptr->glb_selection_lns_lns[bin][index_tau];
	  }
	  
	  if (_index_tt_in_range_(ptr->index_tt_nc_g4, ptr->selection_num_nc, ppt->has_nc_gr4)) {
	    rescaling=ptr->glb_selection_lns_gr4[bin][index_tau];
	  }

	  if (_index_tt_in_range_(ptr->index_tt_nc_g5, ptr->selection_num_nc, ppt->has_nc_gr5)) {
	    rescaling=ptr->glb_selection_lns_gr5[bin][index_tau]*ptr->k[index_md][index_q];
	  }

          /* copy from input array to output array */
          sources[index_tau] *= rescaling;
        }
      }

      /* Intrinsic alignment */
      if ((_index_tt_in_range_(ptr->index_tt_alignment,ptr->selection_num_wl,ppt->has_intrinsic_alignment))) {
        /* bin number associated to particular redshift bin and selection function */
	bin = index_tt - ptr->index_tt_alignment;
	
        /* allocate temporary arrays for storing sources and for calling background */
	class_alloc(pvecback,pba->bg_size*sizeof(double),ptr->error_message);

        /* redefine the time sampling */
	memcpy(tau0_minus_tau,ptr->glb_selection_0_tau0_minus_tau_wl[bin],
	       ptr->glb_selection_0_size_wl[bin]*sizeof(double));

        /* resample the source at those times */
        class_call(transfer_source_resample(ppr,
                                            pba,
                                            ppt,
                                            ptr,
                                            bin,
                                            tau0_minus_tau,
                                            tau_size,
                                            index_md,
                                            tau0,
                                            interpolated_sources,
                                            sources),
                   ptr->error_message,
                   ptr->error_message);

        /* Compute trapezoidal weights for integration over tau */
        class_call(array_trapezoidal_mweights(tau0_minus_tau,
                                              tau_size,
                                              w_trapz,
                                              ptr->error_message),
                   ptr->error_message,
                   ptr->error_message);

        /* compute values of selection function at sampled values of tau */
	selection=ptr->glb_selection_0_wl[bin];
	
	//Compute tracer type
	itr=ptr->selection_tracer_wl[bin];

        /* loop over time and rescale */
        for (index_tau = 0; index_tau < tau_size; index_tau++) {
	  double x=ptr->k[index_md][index_q]*tau0_minus_tau[index_tau];
	  if(x==0)
	    rescaling=0;
	  else {
	    double fred,bI,aI,Dz;

	    /* conformal time */
	    tau = tau0 - tau0_minus_tau[index_tau];
	    
	    /* corresponding background quantities */
	    class_call(background_at_tau(pba,
					 tau,
					 pba->long_info,
					 pba->inter_normal,
					 &last_index,
					 pvecback),
		       pba->error_message,
		       ptr->error_message);
	    
	    /* Clustering bias for density term */
	    z = pba->a_today/pvecback[pba->index_bg_a]-1.;
	    class_call(transfer_interpolate_bias_function(ppr,pba,ppt,ptr,
							  z,itr,
							  ptr->aI_z,ptr->aI_aI,ptr->aI_ddaI,ptr->aI_size,
							  &last_index,&aI),
		       ptr->error_message,ptr->error_message);
	    class_call(transfer_interpolate_bias_function(ppr,pba,ppt,ptr,
							  z,itr,
							  ptr->fred_z,ptr->fred_fred,
							  ptr->fred_ddfred,ptr->fred_size,
							  &last_index,&fred),
		       ptr->error_message,ptr->error_message);
	    
	    /* alignment source =  [- (dz/dtau) W(z)] * b^I * x_red * delta_m(k,tau) /(k*chi)^2 */
	    //This is for the parametrization in 1506.02671
	    Dz=pvecback[pba->index_bg_D]*ppt->inv_growth_0;
	    bI=-aI*(pba->Omega0_b+pba->Omega0_cdm)*0.0134/Dz;
	    //Another option would be just:
	    //bI=aI;

	    rescaling = bI*fred*selection[index_tau]/(x*x);
	  }
	  sources[index_tau] *= rescaling;
        }

        /* deallocate temporary arrays */
	free(pvecback);
      }

      /* lensing potential: eliminate early times, and multiply by selection function */
      if ((_index_tt_in_range_(ptr->index_tt_lensing, ptr->selection_num_wl, ppt->has_lensing_shear))) {
	bin = index_tt - ptr->index_tt_lensing;

        /* redefine the time sampling */
	memcpy(tau0_minus_tau,ptr->glb_selection_lns_tau0_minus_tau_wl[bin],
	       ptr->glb_selection_lns_size_wl[bin]*sizeof(double));

        /* resample the source at those times */
        class_call(transfer_source_resample(ppr,
                                            pba,
                                            ppt,
                                            ptr,
                                            bin,
                                            tau0_minus_tau,
                                            tau_size,
                                            index_md,
                                            tau0,
                                            interpolated_sources,
                                            sources),
                   ptr->error_message,
                   ptr->error_message);

        /* Compute trapezoidal weights for integration over tau */
        class_call(array_trapezoidal_mweights(tau0_minus_tau,
                                              tau_size,
                                              w_trapz,
                                              ptr->error_message),
                   ptr->error_message,
                   ptr->error_message);

        /* loop over time and rescale */
        for (index_tau = 0; index_tau < tau_size; index_tau++) {
          /* lensing source =  - W(tau) (phi(k,tau) + psi(k,tau)) Heaviside(tau-tau_rec)
             with
             psi,phi = metric perturbation in newtonian gauge (phi+psi = Phi_A-Phi_H of Bardeen)
             W = (tau-tau_rec)/(tau_0-tau)/(tau_0-tau_rec)
             H(x) = Heaviside
             (in tau = tau_0, set source = 0 to avoid division by zero;
             regulated anyway by Bessel).
          */
	  rescaling=ptr->glb_selection_lns_wl[bin][index_tau];

          /* copy from input array to output array */
          sources[index_tau] *= rescaling;
        }
      }
    }
  }

  /* case where we do not need to redefine */

  else {

    /* number of sampled time values */
    tau_size = ppt->tau_size;

    /* plain copy from input array to output array */
    memcpy(sources,
           interpolated_sources,
           ppt->tau_size*sizeof(double));

    /* store values of (tau0-tau) */
    for (index_tau=0; index_tau < ppt->tau_size; index_tau++) {
      tau0_minus_tau[index_tau] = tau0 - ppt->tau_sampling[index_tau];
    }

    /* Compute trapezoidal weights for integration over tau */
    class_call(array_trapezoidal_mweights(tau0_minus_tau,
                                          tau_size,
                                          w_trapz,
                                          ptr->error_message),
               ptr->error_message,
               ptr->error_message);
  }

  /* return tau_size value that will be stored in the workspace (the
     workspace wants a double) */

  *tau_size_out = tau_size;

  return _SUCCESS_;

}

/**
 * arbitrarily normalized selection function dN/dz(z,bin)
 *
 * @param ppt                   Input : pointer to perturbation structure
 * @param ptr                   Input : pointer to transfers structure
 * @param bin                   Input : redshift bin number
 * @param z                     Input : one value of redshift
 * @param selection             Output: pointer to selection function
 * @return the error status
 */

int transfer_selection_function_nc(
				   struct precision * ppr,
				   struct perturbs * ppt,
				   struct transfers * ptr,
				   int bin,
				   double z,
				   double * selection) {

  double x;
  double dNdz;
  int last_index;
  short itr=ptr->selection_tracer_nc[bin];

  if((z<ptr->nz_z_nc[itr][0]) || (z>ptr->nz_z_nc[itr][ptr->nz_size_nc[itr]-1])) {
    *selection=0.;
    return _SUCCESS_;
  }
  
  class_call(array_interpolate_spline(
				      ptr->nz_z_nc[itr],
				      ptr->nz_size_nc[itr],
				      ptr->nz_nz_nc[itr],
				      ptr->nz_ddnz_nc[itr],
				      1,
				      z,
				      &last_index,
				      &dNdz,
				      1,
				      ptr->error_message),
	     ptr->error_message,
	     ptr->error_message);
  
  /* difference between z and the bin center (we can take the absolute
     value as long as all selection functions are symmetric around
     x=0) */
  x=fabs(z-ptr->selection_mean_nc[bin]);

  /* gaussian case (the function is anyway normalized later
     automatically, but could not resist to normalize it already
     here) */
  if (ptr->selection_nc[itr]==gaussian) {
    *selection = exp(-0.5*pow(x/ptr->selection_width_nc[bin],2))
      /ptr->selection_width_nc[bin]/sqrt(2.*_PI_);

    *selection *= dNdz;
    return _SUCCESS_;
  }
  
  /* top-hat case, with smoothed edges. The problem with sharp edges
     is that the final result will be affected by random
     noise. Indeed, the values of k at which the transfer functions
     Delta_l(k) are sampled will never coicide with the actual edges
     of the true transfer function (computed with or even without the
     Limber approximation). Hence the integral Cl=\int dk
     Delta_l(k)**2 (...) will be unprecise and will fluctuate randomly
     with the resolution along k. With smooth edges, the problem is
     sloved, and the final Cls become mildly dependent on the
     resolution along k. */

  if (ptr->selection_nc[itr]==tophat) {

    /* selection function, centered on z=mean (i.e. on x=0), equal to
       one around x=0, with tanh step centered on x=width, of width
       delta x = 0.1*width
    */
    if(ptr->has_photoz_nc[itr]==_TRUE_) {
      *selection=erf_mine((ptr->selection_mean_nc[bin]+ptr->selection_width_nc[bin]-z)/(M_SQRT2*ptr->selection_sz_nc[bin]))-
	erf_mine((ptr->selection_mean_nc[bin]-ptr->selection_width_nc[bin]-z)/(M_SQRT2*ptr->selection_sz_nc[bin]));
    }
    else
      *selection=(1.-tanh((x-ptr->selection_width_nc[bin])/(ppr->selection_tophat_edge*ptr->selection_width_nc[bin])))/2.;
    
    *selection *= dNdz;
    return _SUCCESS_;
  }

  /* get here only if selection type was recognized */
  class_stop(ptr->error_message,
             "invalid choice of selection function");

  return _SUCCESS_;
}

/**
 * arbitrarily normalized selection function dN/dz(z,bin)
 *
 * @param ppt                   Input : pointer to perturbation structure
 * @param ptr                   Input : pointer to transfers structure
 * @param bin                   Input : redshift bin number
 * @param z                     Input : one value of redshift
 * @param selection             Output: pointer to selection function
 * @return the error status
 */

int transfer_selection_function_wl(
				   struct precision * ppr,
				   struct perturbs * ppt,
				   struct transfers * ptr,
				   int bin,
				   double z,
				   double * selection) {

  double x;
  double dNdz;
  int last_index;
  short itr=ptr->selection_tracer_wl[bin];

  if((z<ptr->nz_z_wl[itr][0]) || (z>ptr->nz_z_wl[itr][ptr->nz_size_wl[itr]-1])) {
    *selection=0.;
    return _SUCCESS_;
  }
  
  class_call(array_interpolate_spline(
				      ptr->nz_z_wl[itr],
				      ptr->nz_size_wl[itr],
				      ptr->nz_nz_wl[itr],
				      ptr->nz_ddnz_wl[itr],
				      1,
				      z,
				      &last_index,
				      &dNdz,
				      1,
				      ptr->error_message),
	     ptr->error_message,
	     ptr->error_message);
  
  /* difference between z and the bin center (we can take the absolute
     value as long as all selection functions are symmetric around
     x=0) */
  x=fabs(z-ptr->selection_mean_wl[bin]);

  /* gaussian case (the function is anyway normalized later
     automatically, but could not resist to normalize it already
     here) */
  if (ptr->selection_wl[itr]==gaussian) {
    *selection = exp(-0.5*pow(x/ptr->selection_width_wl[bin],2))
      /ptr->selection_width_wl[bin]/sqrt(2.*_PI_);

    *selection *= dNdz;
    return _SUCCESS_;
  }
  
  /* top-hat case, with smoothed edges. The problem with sharp edges
     is that the final result will be affected by random
     noise. Indeed, the values of k at which the transfer functions
     Delta_l(k) are sampled will never coicide with the actual edges
     of the true transfer function (computed with or even without the
     Limber approximation). Hence the integral Cl=\int dk
     Delta_l(k)**2 (...) will be unprecise and will fluctuate randomly
     with the resolution along k. With smooth edges, the problem is
     sloved, and the final Cls become mildly dependent on the
     resolution along k. */

  if (ptr->selection_wl[itr]==tophat) {

    /* selection function, centered on z=mean (i.e. on x=0), equal to
       one around x=0, with tanh step centered on x=width, of width
       delta x = 0.1*width
    */
    if(ptr->has_photoz_wl[itr]==_TRUE_) {
      *selection=erf_mine((ptr->selection_mean_wl[bin]+ptr->selection_width_wl[bin]-z)/(M_SQRT2*ptr->selection_sz_wl[bin]))-
	erf_mine((ptr->selection_mean_wl[bin]-ptr->selection_width_wl[bin]-z)/(M_SQRT2*ptr->selection_sz_wl[bin]));
    }
    else
      *selection=(1.-tanh((x-ptr->selection_width_wl[bin])/(ppr->selection_tophat_edge*ptr->selection_width_wl[bin])))/2.;
    
    *selection *= dNdz;
    return _SUCCESS_;
  }

  /* get here only if selection type was recognized */
  class_stop(ptr->error_message,
             "invalid choice of selection function");

  return _SUCCESS_;
}

/**
 * for sources that need to be mutiplied by a selection function,
 * redefine a finer time sampling in a small range
 *
 * @param ppr                   Input : pointer to precision structure
 * @param pba                   Input : pointer to background structure
 * @param ppt                   Input : pointer to perturbation structure
 * @param ptr                   Input : pointer to transfers structure
 * @param bin                   Input : redshift bin number
 * @param tau0_minus_tau        Output: values of (tau0-tau) at which source are sample
 * @param tau_size              Output: pointer to size of previous array
 * @param index_md            Input : index of mode
 * @param tau0                  Input : time today
 * @param interpolated_sources  Input : interpolated perturbation source
 * @param sources               Output: resampled transfer source
 * @return the error status
 */

int transfer_selection_sampling_nc(
				   struct precision * ppr,
				   struct background * pba,
				   struct perturbs * ppt,
				   struct transfers * ptr,
				   int bin,
				   double * tau0_minus_tau,
				   int tau_size) {
  
  /* running index on time */
  int index_tau;

  /* minimum and maximal value of time in new sampled interval */
  double tau_min,tau_mean,tau_max;

  /* time interval for this bin */
  class_call(transfer_selection_times_nc(ppr,
					 pba,
					 ppt,
					 ptr,
					 bin,
					 &tau_min,
					 &tau_mean,
					 &tau_max),
             ptr->error_message,
             ptr->error_message);

  /* case selection == dirac */
  if (tau_min == tau_max) {
    class_test(tau_size !=1,
               ptr->error_message,
               "for Dirac selection function tau_size should be 1, not %d",tau_size);
    tau0_minus_tau[0] = pba->conformal_age - tau_mean;
  }
  /* for other cases (gaussian, tophat...) define new sampled values
     of (tau0-tau) with even spacing */
  else {
    for (index_tau=0; index_tau<tau_size; index_tau++) {
      tau0_minus_tau[index_tau]=pba->conformal_age-tau_min-((double)index_tau)/((double)tau_size-1.)*(tau_max-tau_min);
    }
  }

  return _SUCCESS_;

}

/**
 * for sources that need to be mutiplied by a selection function,
 * redefine a finer time sampling in a small range
 *
 * @param ppr                   Input : pointer to precision structure
 * @param pba                   Input : pointer to background structure
 * @param ppt                   Input : pointer to perturbation structure
 * @param ptr                   Input : pointer to transfers structure
 * @param bin                   Input : redshift bin number
 * @param tau0_minus_tau        Output: values of (tau0-tau) at which source are sample
 * @param tau_size              Output: pointer to size of previous array
 * @param index_md            Input : index of mode
 * @param tau0                  Input : time today
 * @param interpolated_sources  Input : interpolated perturbation source
 * @param sources               Output: resampled transfer source
 * @return the error status
 */

int transfer_selection_sampling_wl(
				   struct precision * ppr,
				   struct background * pba,
				   struct perturbs * ppt,
				   struct transfers * ptr,
				   int bin,
				   double * tau0_minus_tau,
				   int tau_size) {

  /* running index on time */
  int index_tau;

  /* minimum and maximal value of time in new sampled interval */
  double tau_min,tau_mean,tau_max;

  /* time interval for this bin */
  class_call(transfer_selection_times_wl(ppr,
					 pba,
					 ppt,
					 ptr,
					 bin,
					 &tau_min,
					 &tau_mean,
					 &tau_max),
             ptr->error_message,
             ptr->error_message);

  /* case selection == dirac */
  if (tau_min == tau_max) {
    class_test(tau_size !=1,
               ptr->error_message,
               "for Dirac selection function tau_size should be 1, not %d",tau_size);
    tau0_minus_tau[0] = pba->conformal_age - tau_mean;
  }
  /* for other cases (gaussian, tophat...) define new sampled values
     of (tau0-tau) with even spacing */
  else {
    for (index_tau=0; index_tau<tau_size; index_tau++) {
      tau0_minus_tau[index_tau]=pba->conformal_age-tau_min-((double)index_tau)/((double)tau_size-1.)*(tau_max-tau_min);
    }
  }

  return _SUCCESS_;

}

/**
 * for lensing sources that need to be convolved with a selection
 * function, redefine the sampling within the range extending from the
 * tau_min of the selection function up to tau0
 *
 * @param ppr                   Input : pointer to precision structure
 * @param pba                   Input : pointer to background structure
 * @param ppt                   Input : pointer to perturbation structure
 * @param ptr                   Input : pointer to transfers structure
 * @param bin                   Input : redshift bin number
 * @param tau0_minus_tau        Output: values of (tau0-tau) at which source are sample
 * @param tau_size              Output: pointer to size of previous array
 * @param index_md            Input : index of mode
 * @param tau0                  Input : time today
 * @param interpolated_sources  Input : interpolated perturbation source
 * @param sources               Output: resampled transfer source
 * @return the error status
 */

int transfer_lensing_sampling_nc(
				 struct precision * ppr,
				 struct background * pba,
				 struct perturbs * ppt,
				 struct transfers * ptr,
				 int bin,
				 double tau0,
				 double * tau0_minus_tau,
				 int tau_size) {

  /* running index on time */
  int index_tau;

  /* minimum and maximal value of time in new sampled interval */
  double tau_min,tau_mean,tau_max;

  /* time interval for this bin */
  class_call(transfer_selection_times_nc(ppr,
					 pba,
					 ppt,
					 ptr,
					 bin,
					 &tau_min,
					 &tau_mean,
					 &tau_max),
             ptr->error_message,
             ptr->error_message);

  for (index_tau=0; index_tau<tau_size; index_tau++) {
    //tau0_minus_tau[index_tau]=pba->conformal_age-tau_min-((double)index_tau)/((double)tau_size-1.)*(tau0-tau_min);
    tau0_minus_tau[index_tau]=((double)(tau_size-1-index_tau))/((double)(tau_size-1))*(tau0-tau_min);
  }

  return _SUCCESS_;

}

/**
 * for lensing sources that need to be convolved with a selection
 * function, redefine the sampling within the range extending from the
 * tau_min of the selection function up to tau0
 *
 * @param ppr                   Input : pointer to precision structure
 * @param pba                   Input : pointer to background structure
 * @param ppt                   Input : pointer to perturbation structure
 * @param ptr                   Input : pointer to transfers structure
 * @param bin                   Input : redshift bin number
 * @param tau0_minus_tau        Output: values of (tau0-tau) at which source are sample
 * @param tau_size              Output: pointer to size of previous array
 * @param index_md            Input : index of mode
 * @param tau0                  Input : time today
 * @param interpolated_sources  Input : interpolated perturbation source
 * @param sources               Output: resampled transfer source
 * @return the error status
 */

int transfer_lensing_sampling_wl(
				 struct precision * ppr,
				 struct background * pba,
				 struct perturbs * ppt,
				 struct transfers * ptr,
				 int bin,
				 double tau0,
				 double * tau0_minus_tau,
				 int tau_size) {
  
  /* running index on time */
  int index_tau;

  /* minimum and maximal value of time in new sampled interval */
  double tau_min,tau_mean,tau_max;

  /* time interval for this bin */
  class_call(transfer_selection_times_wl(ppr,
					 pba,
					 ppt,
					 ptr,
					 bin,
					 &tau_min,
					 &tau_mean,
					 &tau_max),
             ptr->error_message,
             ptr->error_message);

  for (index_tau=0; index_tau<tau_size; index_tau++) {
    //tau0_minus_tau[index_tau]=pba->conformal_age-tau_min-((double)index_tau)/((double)tau_size-1.)*(tau0-tau_min);
    tau0_minus_tau[index_tau]=((double)(tau_size-1-index_tau))/((double)(tau_size-1))*(tau0-tau_min);
  }

  return _SUCCESS_;

}


/**
 * for sources that need to be mutiplied by a selection function,
 * redefine a finer time sampling in a small range, and resample the
 * perturbation sources at the new value by linear interpolation
 *
 * @param ppr                   Input : pointer to precision structure
 * @param pba                   Input : pointer to background structure
 * @param ppt                   Input : pointer to perturbation structure
 * @param ptr                   Input : pointer to transfers structure
 * @param bin                   Input : redshift bin number
 * @param tau0_minus_tau        Output: values of (tau0-tau) at which source are sample
 * @param tau_size              Output: pointer to size of previous array
 * @param index_md            Input : index of mode
 * @param tau0                  Input : time today
 * @param interpolated_sources  Input : interpolated perturbation source
 * @param sources               Output: resampled transfer source
 * @return the error status
 */

int transfer_source_resample(
                             struct precision * ppr,
                             struct background * pba,
                             struct perturbs * ppt,
                             struct transfers * ptr,
                             int bin,
                             double * tau0_minus_tau,
                             int tau_size,
                             int index_md,
                             double tau0,
                             double * interpolated_sources,
                             double * sources) {

  /* running index on time */
  int index_tau;

  /* array of values of source */
  double * source_at_tau;

  /* array of source values for a given time and for all k's */
  class_alloc(source_at_tau,
              sizeof(double),
              ptr->error_message);

  /* interpolate the sources linearily at the new time values */
  for (index_tau=0; index_tau<tau_size; index_tau++) {

    class_call(array_interpolate_two(ppt->tau_sampling,
                                     1,
                                     0,
                                     interpolated_sources,
                                     1,
                                     ppt->tau_size,
                                     tau0-tau0_minus_tau[index_tau],
                                     source_at_tau,
                                     1,
                                     ptr->error_message),
               ptr->error_message,
               ptr->error_message);

    /* copy the new values in the output sources array */
    sources[index_tau] = source_at_tau[0];
  }

  /* deallocate the temporary array */
  free(source_at_tau);

  return _SUCCESS_;

}

/**
 * for each selection function, compute the min, mean and max values
 * of conformal time (associated to the min, mean and max values of
 * redshift specified by the user)
 *
 * @param ppr                   Input : pointer to precision structure
 * @param pba                   Input : pointer to background structure
 * @param ppt                   Input : pointer to perturbation structure
 * @param ptr                   Input : pointer to transfers structure
 * @param bin                   Input : redshift bin number
 * @param tau_min               Output: smallest time in the selection interval
 * @param tau_mean              Output: time corresponding to z_mean
 * @param tau_max               Output: largest time in the selection interval
 * @return the error status
 */

int transfer_selection_times_nc(
				struct precision * ppr,
				struct background * pba,
				struct perturbs * ppt,
				struct transfers * ptr,
				int bin,
				double * tau_min,
				double * tau_mean,
				double * tau_max) {

  /* a value of redshift */
  double z=0.;
  short itr=ptr->selection_tracer_nc[bin];

  /* lower edge of time interval for this bin */

  if (ptr->selection_nc[itr]==gaussian) {
    z = ptr->selection_mean_nc[bin]+ptr->selection_width_nc[bin]*ppr->selection_cut_at_sigma;
  }
  if (ptr->selection_nc[itr]==tophat) {
    if(ptr->has_photoz_nc[itr]==_TRUE_)
      z = ptr->selection_mean_nc[bin]+ptr->selection_width_nc[bin]+3*ptr->selection_sz_nc[bin];
    else 
      z = ptr->selection_mean_nc[bin]+(1.+ppr->selection_cut_at_sigma*ppr->selection_tophat_edge)*ptr->selection_width_nc[bin];
  }

  class_call(background_tau_of_z(pba,
                                 z,
                                 tau_min),
             pba->error_message,
             ppt->error_message);

  /* higher edge of time interval for this bin */

  if (ptr->selection_nc[itr]==gaussian) {
    z = MAX(ptr->selection_mean_nc[bin]-ptr->selection_width_nc[bin]*ppr->selection_cut_at_sigma,0.);
  }
  if (ptr->selection_nc[itr]==tophat) {
    if(ptr->has_photoz_nc[itr]==_TRUE_)
      z = MAX(ptr->selection_mean_nc[bin]-ptr->selection_width_nc[bin]-3*ptr->selection_sz_nc[bin],0.);
    else
      z = MAX(ptr->selection_mean_nc[bin]-(1.+ppr->selection_cut_at_sigma*ppr->selection_tophat_edge)*ptr->selection_width_nc[bin],0.);
  }

  class_call(background_tau_of_z(pba,
                                 z,
                                 tau_max),
             pba->error_message,
             ppt->error_message);

  /* central value of time interval for this bin */

  z = MAX(ptr->selection_mean_nc[bin],0.);

  class_call(background_tau_of_z(pba,
                                 z,
                                 tau_mean),
             pba->error_message,
             ppt->error_message);

  return _SUCCESS_;

}

/**
 * for each selection function, compute the min, mean and max values
 * of conformal time (associated to the min, mean and max values of
 * redshift specified by the user)
 *
 * @param ppr                   Input : pointer to precision structure
 * @param pba                   Input : pointer to background structure
 * @param ppt                   Input : pointer to perturbation structure
 * @param ptr                   Input : pointer to transfers structure
 * @param bin                   Input : redshift bin number
 * @param tau_min               Output: smallest time in the selection interval
 * @param tau_mean              Output: time corresponding to z_mean
 * @param tau_max               Output: largest time in the selection interval
 * @return the error status
 */

int transfer_selection_times_wl(
				struct precision * ppr,
				struct background * pba,
				struct perturbs * ppt,
				struct transfers * ptr,
				int bin,
				double * tau_min,
				double * tau_mean,
				double * tau_max) {
  
  /* a value of redshift */
  double z=0.;
  short itr=ptr->selection_tracer_wl[bin];

  /* lower edge of time interval for this bin */

  if (ptr->selection_wl[itr]==gaussian) {
    z = ptr->selection_mean_wl[bin]+ptr->selection_width_wl[bin]*ppr->selection_cut_at_sigma;
  }
  if (ptr->selection_wl[itr]==tophat) {
    if(ptr->has_photoz_wl[itr]==_TRUE_)
      z = ptr->selection_mean_wl[bin]+ptr->selection_width_wl[bin]+3*ptr->selection_sz_wl[bin];
    else 
      z = ptr->selection_mean_wl[bin]+(1.+ppr->selection_cut_at_sigma*ppr->selection_tophat_edge)*ptr->selection_width_wl[bin];
  }

  class_call(background_tau_of_z(pba,
                                 z,
                                 tau_min),
             pba->error_message,
             ppt->error_message);

  /* higher edge of time interval for this bin */

  if (ptr->selection_wl[itr]==gaussian) {
    z = MAX(ptr->selection_mean_wl[bin]-ptr->selection_width_wl[bin]*ppr->selection_cut_at_sigma,0.);
  }
  if (ptr->selection_wl[itr]==tophat) {
    if(ptr->has_photoz_wl[itr]==_TRUE_)
      z = MAX(ptr->selection_mean_wl[bin]-ptr->selection_width_wl[bin]-3*ptr->selection_sz_wl[bin],0.);
    else
      z = MAX(ptr->selection_mean_wl[bin]-(1.+ppr->selection_cut_at_sigma*ppr->selection_tophat_edge)*ptr->selection_width_wl[bin],0.);
  }

  class_call(background_tau_of_z(pba,
                                 z,
                                 tau_max),
             pba->error_message,
             ppt->error_message);

  /* central value of time interval for this bin */

  z = MAX(ptr->selection_mean_wl[bin],0.);

  class_call(background_tau_of_z(pba,
                                 z,
                                 tau_mean),
             pba->error_message,
             ppt->error_message);

  return _SUCCESS_;

}

/**
 * compute and normalise selection function for a set of time values
 *
 *
 * @param pba                   Input : pointer to background structure
 * @param ppt                   Input : pointer to perturbation structure
 * @param ptr                   Input : pointer to transfers structure
 * @param selection             Output: normalized selection function
 * @param tau0_minus_tau        Input : values of (tau0-tau) at which source are sample
 * @param w_trapz               Input : trapezoidal weights for integration over tau
 * @param tau_size              Input : size of previous two arrays
 * @param pvecback              Input : allocated array of background values
 * @param tau_0                 Input : time today
 * @param bin                   Input : redshift bin number
 * @return the error status
 */

int transfer_selection_compute_nc(
				  struct precision * ppr,
				  struct background * pba,
				  struct perturbs * ppt,
				  struct transfers * ptr,
				  double * selection,
				  double * tau0_minus_tau,
				  double * w_trapz,
				  int tau_size,
				  double * pvecback,
				  double tau0,
				  int bin) {

  /* running index over time */
  int index_tau;

  /* running value of time */
  double tau;

  /* used for normalizing the selection to one */
  double norm;

  /* used for calling background_at_tau() */
  int last_index;

  /* runnign value of redshift */
  double z;

  if (tau_size > 1) {

    /* loop over time */
    for (index_tau = 0; index_tau < tau_size; index_tau++) {
      if(tau0_minus_tau[index_tau]<1E-5)
	tau0_minus_tau[index_tau]=0;
      /* running value of time */
      tau = tau0 - tau0_minus_tau[index_tau];

      /* get background quantitites at this time */
      class_call(background_at_tau(pba,
                                   tau,
                                   pba->long_info,
                                   pba->inter_normal,
                                   &last_index,
                                   pvecback),
                 pba->error_message,
                 ptr->error_message);

      /* infer redhsift */
      z = pba->a_today/pvecback[pba->index_bg_a]-1.;

      /* get corresponding dN/dz(z,bin) */
      class_call(transfer_selection_function_nc(ppr,
						ppt,
						ptr,
						bin,
						z,
						&(selection[index_tau])),
                 ptr->error_message,
                 ptr->error_message);

      /* get corresponding dN/dtau = dN/dz * dz/dtau = dN/dz * H */
      selection[index_tau] *= pvecback[pba->index_bg_H];

    }

    /* compute norm = \int W(tau) dtau */
    class_call(array_trapezoidal_integral(selection,
                                          tau_size,
                                          w_trapz,
                                          &norm,
                                          ptr->error_message),
               ptr->error_message,
               ptr->error_message);


    /* divide W by norm so that \int W(tau) dtau = 1 */
    for (index_tau = 0; index_tau < tau_size; index_tau++) {
      selection[index_tau]/=norm;
    }

  }

  /* trivial case: dirac distribution */
  else {
    selection[0] = 1.;
  }

  return _SUCCESS_;
}

/**
 * compute and normalise selection function for a set of time values
 *
 *
 * @param pba                   Input : pointer to background structure
 * @param ppt                   Input : pointer to perturbation structure
 * @param ptr                   Input : pointer to transfers structure
 * @param selection             Output: normalized selection function
 * @param tau0_minus_tau        Input : values of (tau0-tau) at which source are sample
 * @param w_trapz               Input : trapezoidal weights for integration over tau
 * @param tau_size              Input : size of previous two arrays
 * @param pvecback              Input : allocated array of background values
 * @param tau_0                 Input : time today
 * @param bin                   Input : redshift bin number
 * @return the error status
 */

int transfer_selection_compute_wl(
				  struct precision * ppr,
				  struct background * pba,
				  struct perturbs * ppt,
				  struct transfers * ptr,
				  double * selection,
				  double * tau0_minus_tau,
				  double * w_trapz,
				  int tau_size,
				  double * pvecback,
				  double tau0,
				  int bin) {
  
  /* running index over time */
  int index_tau;

  /* running value of time */
  double tau;

  /* used for normalizing the selection to one */
  double norm;

  /* used for calling background_at_tau() */
  int last_index;

  /* runnign value of redshift */
  double z;

  if (tau_size > 1) {

    /* loop over time */
    for (index_tau = 0; index_tau < tau_size; index_tau++) {
      if(tau0_minus_tau[index_tau]<1E-5)
	tau0_minus_tau[index_tau]=0;
      /* running value of time */
      tau = tau0 - tau0_minus_tau[index_tau];

      /* get background quantitites at this time */
      class_call(background_at_tau(pba,
                                   tau,
                                   pba->long_info,
                                   pba->inter_normal,
                                   &last_index,
                                   pvecback),
                 pba->error_message,
                 ptr->error_message);

      /* infer redhsift */
      z = pba->a_today/pvecback[pba->index_bg_a]-1.;

      /* get corresponding dN/dz(z,bin) */
      class_call(transfer_selection_function_wl(ppr,
						ppt,
						ptr,
						bin,
						z,
						&(selection[index_tau])),
                 ptr->error_message,
                 ptr->error_message);

      /* get corresponding dN/dtau = dN/dz * dz/dtau = dN/dz * H */
      selection[index_tau] *= pvecback[pba->index_bg_H];

    }

    /* compute norm = \int W(tau) dtau */
    class_call(array_trapezoidal_integral(selection,
                                          tau_size,
                                          w_trapz,
                                          &norm,
                                          ptr->error_message),
               ptr->error_message,
               ptr->error_message);


    /* divide W by norm so that \int W(tau) dtau = 1 */
    for (index_tau = 0; index_tau < tau_size; index_tau++) {
      selection[index_tau]/=norm;
    }

  }

  /* trivial case: dirac distribution */
  else {
    selection[0] = 1.;
  }

  return _SUCCESS_;
}

int transfer_interpolate_bias_function(
				       struct precision * ppr,
				       struct background * pba,
				       struct perturbs * ppt,
				       struct transfers * ptr,
				       double z,
				       int itr,
				       double **z_arr,
				       double **fz_arr,
				       double **ddfz_arr,
				       int *fz_size,
				       int *last_index,
				       double *bias_out)
{
  double bias;
  class_test((z<z_arr[itr][0]) || (z>z_arr[itr][fz_size[itr]-1]),
	     ptr->error_message,
	     "One of your input bias functions only covers the redshift"
	     " range [%f : %f]. However I need z=%f",
	     z_arr[itr][0],z_arr[itr][fz_size[itr]-1],z);
  
  class_call(array_interpolate_spline(z_arr[itr],
				      fz_size[itr],
				      fz_arr[itr],
				      ddfz_arr[itr],
				      1,
				      z,
				      last_index,
				      &bias,
				      1,
				      ptr->error_message),
	     ptr->error_message,ptr->error_message);

  *bias_out=bias;

  return _SUCCESS_;
}

/**
 * This routine computes the transfer functions \f$ \Delta_l^{X} (k) \f$)
 * as a function of wavenumber k for a given mode, initial condition,
 * type and multipole l passed in input.
 *
 * For a given value of k, the transfer function is infered from
 * the source function (passed in input in the array interpolated_sources)
 * and from Bessel functions (passed in input in the bessels structure),
 * either by convolving them along tau, or by a Limber appoximation.
 * This elementary task is distributed either to transfer_integrate()
 * or to transfer_limber(). The task of this routine is mainly to
 * loop over k values, and to decide at which k_max the calculation can
 * be stopped, according to some approximation scheme designed to find a
 * compromise between execution time and precision. The approximation scheme
 * is defined by parameters in bthe precision structure.
 *
 * @param ppr                   Input : pointer to precision structure
 * @param ppt                   Input : pointer to perturbation structure
 * @param ptr                   Input/output : pointer to transfers structure (result stored there)
 * @param tau0                  Input : conformal time today
 * @param tau_rec               Input : conformal time at recombination
 * @param index_md            Input : index of mode
 * @param index_ic              Input : index of initial condition
 * @param index_tt              Input : index of type of transfer
 * @param index_l               Input : index of multipole
 * @param interpolated_sources  Input : array containing the sources
 * @param ptw                   Input : pointer to transfer_workspace structure (allocated in transfer_init() to avoid numerous reallocation)
 * @return the error status
 */

int transfer_compute_for_each_l(
                                struct transfer_workspace * ptw,
                                struct precision * ppr,
                                struct perturbs * ppt,
                                struct transfers * ptr,
                                int index_q,
                                int index_md,
                                int index_ic,
                                int index_tt,
                                int index_l,
                                double l,
                                double q_max_bessel,
                                radial_function_type radial_type
                                ){

  /** Summary: */

  /** - define local variables */

  /* current wavenumber value */
  double q,k;

  /* value of transfer function */
  double transfer_function;

  /* whether to use the Limber approximation */
  short use_limber;

  /* return zero tranbsfer function if l is above l_max */
  if (index_l >= ptr->l_size_tt[index_md][index_tt]) {

    ptr->transfer[index_md][((index_ic * ptr->tt_size[index_md] + index_tt)
                             * ptr->l_size[index_md] + index_l)
                            * ptr->q_size + index_q] = 0.;
    return _SUCCESS_;
  }

  q = ptr->q[index_q];
  k = ptr->k[index_md][index_q];

  if (ptr->transfer_verbose > 3)
    printf("Compute transfer for l=%d type=%d\n",(int)l,index_tt);

  class_call(transfer_use_limber(ppr,
                                 ppt,
                                 ptr,
                                 q_max_bessel,
                                 index_md,
                                 index_tt,
                                 q,
                                 l,
                                 &use_limber),
             ptr->error_message,
             ptr->error_message);

  if (use_limber == _TRUE_) {

    class_call(transfer_limber(ptr,
                               ptw,
                               index_md,
                               index_q,
                               l,
                               q,
                               radial_type,
                               &transfer_function),
               ptr->error_message,
               ptr->error_message);

  }
  else {
    class_call(transfer_integrate(
                                  ppt,
                                  ptr,
                                  ptw,
                                  index_q,
                                  index_md,
                                  index_tt,
                                  l,
                                  index_l,
                                  k,
                                  radial_type,
                                  &transfer_function
                                  ),
               ptr->error_message,
               ptr->error_message);
  }

  /* store transfer function in transfer structure */
  ptr->transfer[index_md][((index_ic * ptr->tt_size[index_md] + index_tt)
                           * ptr->l_size[index_md] + index_l)
                          * ptr->q_size + index_q]
    = transfer_function;

  return _SUCCESS_;

}

int transfer_use_limber(
                        struct precision * ppr,
                        struct perturbs * ppt,
                        struct transfers * ptr,
                        double q_max_bessel,
                        int index_md,
                        int index_tt,
                        double q,
                        double l,
                        short * use_limber) {


  /* criterium for chosing between integration and Limber
     must be implemented here */

  *use_limber = _FALSE_;

  if (q>q_max_bessel) {
    *use_limber = _TRUE_;
  }
  else {

    if (_scalars_) {

      //TBC: in principle the Limber condition should be adapted to account for curvature effects

      if ((ppt->has_cl_cmb_lensing_potential == _TRUE_) && (index_tt == ptr->index_tt_lcmb) && (l>ppr->l_switch_limber)) {
        *use_limber = _TRUE_;
      }
      if (_index_tt_in_range_(ptr->index_tt_density, ptr->selection_num_nc, ppt->has_nc_density) && (l>=ppr->l_switch_limber_for_cl_density)) {
        *use_limber = _TRUE_;
      }
      if (_index_tt_in_range_(ptr->index_tt_rsd,     ptr->selection_num_nc, ppt->has_nc_rsd1) && (l>=ppr->l_switch_limber_for_cl_density)) {
        *use_limber = _TRUE_;
      }
      if (_index_tt_in_range_(ptr->index_tt_d0,      ptr->selection_num_nc, ppt->has_nc_rsd2) && (l>=ppr->l_switch_limber_for_cl_density)) {
        *use_limber = _TRUE_;
      }
      if (_index_tt_in_range_(ptr->index_tt_d1,      ptr->selection_num_nc, ppt->has_nc_rsd3) && (l>=ppr->l_switch_limber_for_cl_density)) {
        *use_limber = _TRUE_;
      }
      if (_index_tt_in_range_(ptr->index_tt_nc_lens, ptr->selection_num_nc, ppt->has_nc_lens) && (l>=ppr->l_switch_limber_for_cl_density)) {
        *use_limber = _TRUE_;
      }
      if (_index_tt_in_range_(ptr->index_tt_nc_g1, ptr->selection_num_nc, ppt->has_nc_gr1) && (l>=ppr->l_switch_limber_for_cl_density)) {
        *use_limber = _TRUE_;
      }
      if (_index_tt_in_range_(ptr->index_tt_nc_g2, ptr->selection_num_nc, ppt->has_nc_gr2) && (l>=ppr->l_switch_limber_for_cl_density)) {
        *use_limber = _TRUE_;
      }
      if (_index_tt_in_range_(ptr->index_tt_nc_g3, ptr->selection_num_nc, ppt->has_nc_gr3) && (l>=ppr->l_switch_limber_for_cl_density)) {
        *use_limber = _TRUE_;
      }
      if (_index_tt_in_range_(ptr->index_tt_nc_g4, ptr->selection_num_nc, ppt->has_nc_gr4) && (l>=ppr->l_switch_limber_for_cl_density)) {
        *use_limber = _TRUE_;
      }
      if (_index_tt_in_range_(ptr->index_tt_nc_g5, ptr->selection_num_nc, ppt->has_nc_gr5) && (l>=ppr->l_switch_limber_for_cl_density)) {
        *use_limber = _TRUE_;
      }
      if (_index_tt_in_range_(ptr->index_tt_lensing, ptr->selection_num_wl, ppt->has_lensing_shear) && (l>=ppr->l_switch_limber_for_cl_lensing)) {
        *use_limber = _TRUE_;
      }
      if (_index_tt_in_range_(ptr->index_tt_alignment, ptr->selection_num_wl, ppt->has_intrinsic_alignment) && (l>=ppr->l_switch_limber_for_cl_lensing)) {
        *use_limber = _TRUE_;
      }
    }
  }

  return _SUCCESS_;
}

/**
 * This routine computes the transfer functions \f$ \Delta_l^{X} (k) \f$)
 * for each mode, initial condition, type, multipole l and wavenumber k,
 * by convolving  the source function (passed in input in the array
 * interpolated_sources) with Bessel functions (passed in input in the
 * bessels structure).
 *
 * @param ppt                   Input : pointer to perturbation structure
 * @param ptr                   Input : pointer to transfers structure
 * @param tau0                  Input : conformal time today
 * @param tau_rec               Input : conformal time at recombination
 * @param index_md            Input : index of mode
 * @param index_tt              Input : index of type
 * @param index_l               Input : index of multipole
 * @param index_q               Input : index of wavenumber
 * @param interpolated_sources  Input: array of interpolated sources
 * @param ptw                   Input : pointer to transfer_workspace structure (allocated in transfer_init() to avoid numerous reallocation)
 * @param trsf                  Output: transfer function \f$ \Delta_l(k) \f$
 * @return the error status
 */

int transfer_integrate(
                       struct perturbs * ppt,
                       struct transfers * ptr,
                       struct transfer_workspace *ptw,
                       int index_q,
                       int index_md,
                       int index_tt,
                       double l,
                       int index_l,
                       double k,
                       radial_function_type radial_type,
                       double * trsf
                       ) {

  /** Summary: */

  /** - define local variables */

  double * tau0_minus_tau = ptw->tau0_minus_tau;
  double * w_trapz = ptw->w_trapz;
  double * sources = ptw->sources;

  /* minimum value of \f$ (\tau0-\tau) \f$ at which \f$ j_l(k[\tau_0-\tau]) \f$ is known, given that \f$ j_l(x) \f$ is sampled above some finite value \f$ x_{\min} \f$ (below which it can be approximated by zero) */
  double tau0_minus_tau_min_bessel;

  /* index in the source's tau list corresponding to the last point in the overlapping region between sources and bessels. Also the index of possible Bessel truncation. */
  int index_tau_max, index_tau_max_Bessel;

  double bessel, *radial_function;

  double x_turning_point;

  /** - find minimum value of (tau0-tau) at which \f$ j_l(k[\tau_0-\tau]) \f$ is known, given that \f$ j_l(x) \f$ is sampled above some finite value \f$ x_{\min} \f$ (below which it can be approximated by zero) */
  //tau0_minus_tau_min_bessel = x_min_l/k; /* segmentation fault impossible, checked before that k != 0 */
  if (ptw->sgnK==0){
    tau0_minus_tau_min_bessel = ptw->pBIS->chi_at_phimin[index_l]/k; /* segmentation fault impossible, checked before that k != 0 */
  }
  else{

    if (index_q < ptr->index_q_flat_approximation) {

      tau0_minus_tau_min_bessel = ptw->HIS.chi_at_phimin[index_l]/sqrt(ptw->sgnK*ptw->K);

    }
    else {

      tau0_minus_tau_min_bessel = ptw->pBIS->chi_at_phimin[index_l]/sqrt(ptw->sgnK*ptw->K);

      if (ptw->sgnK == 1) {
        x_turning_point = asin(sqrt(l*(l+1.))/ptr->q[index_q]*sqrt(ptw->sgnK*ptw->K));
        tau0_minus_tau_min_bessel *= x_turning_point/sqrt(l*(l+1.));
      }
      else {
        x_turning_point = asinh(sqrt(l*(l+1.))/ptr->q[index_q]*sqrt(ptw->sgnK*ptw->K));
        tau0_minus_tau_min_bessel *= x_turning_point/sqrt(l*(l+1.));
      }
    }
  }
  /** - if there is no overlap between the region in which bessels and sources are non-zero, return zero */
  if (tau0_minus_tau_min_bessel >= tau0_minus_tau[0]) {
    *trsf = 0.;
    return _SUCCESS_;
  }

  /** - if there is an overlap: */

  /** -> trivial case: the source is a Dirac function and is sampled in only one point */
  if (ptw->tau_size == 1) {

    class_call(transfer_radial_function(
                                        ptw,
                                        ppt,
                                        ptr,
                                        k,
                                        index_q,
                                        index_l,
                                        1,
                                        &bessel,
                                        radial_type
                                        ),
               ptr->error_message,
               ptr->error_message);

    *trsf = sources[0] * bessel;
    return _SUCCESS_;
  }

  /** -> other cases */

  /** (a) find index in the source's tau list corresponding to the last point in the overlapping region. After this step, index_tau_max can be as small as zero, but not negative. */
  index_tau_max = ptw->tau_size-1;
  while (tau0_minus_tau[index_tau_max] < tau0_minus_tau_min_bessel)
    index_tau_max--;
  /* Set index so we know if the truncation of the convolution integral is due to Bessel and not
     due to the source. */
  index_tau_max_Bessel = index_tau_max;

  /** (b) the source function can vanish at large $\f \tau \f$. Check if further points can be eliminated. After this step and if we did not return a null transfer function, index_tau_max can be as small as zero, but not negative. */
  while (sources[index_tau_max] == 0.) {
    index_tau_max--;
    if (index_tau_max < 0) {
      *trsf = 0.;
      return _SUCCESS_;
    }
  }

  if (ptw->neglect_late_source == _TRUE_) {

    while (tau0_minus_tau[index_tau_max] < ptw->tau0_minus_tau_cut) {
      index_tau_max--;
      if (index_tau_max < 0) {
        *trsf = 0.;
        return _SUCCESS_;
      }
    }
  }

  /** Compute the radial function: */
  class_alloc(radial_function,sizeof(double)*(index_tau_max+1),ptr->error_message);

  class_call(transfer_radial_function(
                                      ptw,
                                      ppt,
                                      ptr,
                                      k,
                                      index_q,
                                      index_l,
                                      index_tau_max+1,
                                      radial_function,
                                      radial_type
                                      ),
             ptr->error_message,
             ptr->error_message);

  /** Now we do most of the convolution integral: */
  class_call(array_trapezoidal_convolution(sources,
                                           radial_function,
                                           index_tau_max+1,
                                           w_trapz,
                                           trsf,
                                           ptr->error_message),
             ptr->error_message,
             ptr->error_message);

  /** This integral is correct for the case where no truncation has
      occured. If it has been truncated at some index_tau_max because
      f[index_tau_max+1]==0, it is still correct. The 'mistake' in using
      the wrong weight w_trapz[index_tau_max] is exactly compensated by the
      triangle we miss. However, for the Bessel cut off, we must subtract the
      wrong triangle and add the correct triangle */
  if ((index_tau_max!=(ptw->tau_size-1))&&(index_tau_max==index_tau_max_Bessel)){
    //Bessel truncation
    *trsf -= 0.5*(tau0_minus_tau[index_tau_max+1]-tau0_minus_tau_min_bessel)*
      radial_function[index_tau_max]*sources[index_tau_max];
  }


  free(radial_function);
  return _SUCCESS_;
}

/**
 * This routine computes the transfer functions \f$ \Delta_l^{X} (k) \f$)
 * for each mode, initial condition, type, multipole l and wavenumber k,
 * by using the Limber approximation, i.e by evaluating the source function
 * (passed in input in the array interpolated_sources) at a single value of
 * tau (the Bessel function being approximated as a Dirac distribution)
 *
 * @param ptr                   Input : pointer to transfers structure
 * @param ptw                   Input : pointer to transfer workspace structure
 * @param index_md              Input : index of mode
 * @param index_l               Input : index of multipole
 * @param index_q               Input : index of wavenumber
 * @param radial_type           Input : type of radial (Bessel) functions to convolve with
 * @param trsf                  Output: transfer function \f$ \Delta_l(k) \f$
 * @return the error status
 */

int transfer_limber(
                    struct transfers * ptr,
                    struct transfer_workspace * ptw,
                    int index_md,
                    int index_q,
                    double l,
                    double q,
                    radial_function_type radial_type,
                    double * trsf
                    ){

  /** Summary: */

  /** - define local variables */

  /* interpolated source and its derivatives at this value */
  double S, Sp, Sm;

  double x_limber=0.;
  double tau0_minus_tau_limber=0.;
  double IPhiFlat = 0.;

  if (radial_type == SCALAR_TEMPERATURE_0) {

    /** - get k, l and infer tau such that k(tau0-tau)=l+1/2;
        check that tau is in appropriate range */

    if (ptw->sgnK == 0) {
      tau0_minus_tau_limber = (l+0.5)/q;
    }
    else if (ptw->sgnK == 1) {
      x_limber = asin(sqrt(l*(l+1.))/q*sqrt(ptw->K));
      tau0_minus_tau_limber = x_limber/sqrt(ptw->K);
    }
    else if (ptw->sgnK == -1) {
      x_limber = asinh((l+0.5)/q*sqrt(-ptw->K));
      tau0_minus_tau_limber = x_limber/sqrt(-ptw->K);
    }

    if ((tau0_minus_tau_limber > ptw->tau0_minus_tau[0]) ||
        (tau0_minus_tau_limber < ptw->tau0_minus_tau[ptw->tau_size-1])) {
      *trsf = 0.;
      return _SUCCESS_;
    }

    class_call(transfer_limber_interpolate(ptr,
                                           ptw->tau0_minus_tau,
                                           ptw->sources,
                                           ptw->tau_size,
                                           tau0_minus_tau_limber,
                                           &S),
               ptr->error_message,
               ptr->error_message);

    /** - get transfer = source * sqrt(pi/(2l+1))/q
        = source*[tau0-tau] * sqrt(pi/(2l+1))/(l+1/2)
    */

    IPhiFlat = sqrt(_PI_/(2.*l))*(1.-0.25/l+1./32./(l*l));

    *trsf = IPhiFlat*S;

    if (ptw->sgnK == 0) {
      *trsf /= (l+0.5);
    }
    else {
      *trsf *= pow(1.-ptw->K*l*l/q/q,-1./4.)/(tau0_minus_tau_limber*q);
    }

  }

  else if (radial_type == SCALAR_TEMPERATURE_1) {

    if (((l+1.5)/q > ptw->tau0_minus_tau[0]) ||
        ((l-0.5)/q < ptw->tau0_minus_tau[ptw->tau_size-1])) {
      *trsf = 0.;
      return _SUCCESS_;
    }

    class_call(transfer_limber_interpolate(ptr,
                                           ptw->tau0_minus_tau,
                                           ptw->sources,
                                           ptw->tau_size,
                                           (l+1.5)/q,
                                           &Sp),
               ptr->error_message,
               ptr->error_message);

    class_call(transfer_limber_interpolate(ptr,
                                           ptw->tau0_minus_tau,
                                           ptw->sources,
                                           ptw->tau_size,
                                           (l-0.5)/q,
                                           &Sm),
               ptr->error_message,
               ptr->error_message);

    *trsf =
      -sqrt(_PI_/(2.*l+3.))*Sp/(l+1.5) * (l+1.)/(2.*l+1)
      +sqrt(_PI_/(2.*l-1.))*Sm/(l-0.5) * l/(2.*l+1.);

  }

  else if (radial_type == NC_RSD) {

    if (((l+2.5)/q > ptw->tau0_minus_tau[0]) ||
        ((l-1.5)/q < ptw->tau0_minus_tau[ptw->tau_size-1])) {
      *trsf = 0.;
      return _SUCCESS_;
    }

    class_call(transfer_limber_interpolate(ptr,
                                           ptw->tau0_minus_tau,
                                           ptw->sources,
                                           ptw->tau_size,
                                           (l+2.5)/q,
                                           &Sp),
               ptr->error_message,
               ptr->error_message);

    class_call(transfer_limber_interpolate(ptr,
                                           ptw->tau0_minus_tau,
                                           ptw->sources,
                                           ptw->tau_size,
                                           (l-1.5)/q,
                                           &Sm),
               ptr->error_message,
               ptr->error_message);

    class_call(transfer_limber_interpolate(ptr,
                                           ptw->tau0_minus_tau,
                                           ptw->sources,
                                           ptw->tau_size,
                                           (l+0.5)/q,
                                           &S),
               ptr->error_message,
               ptr->error_message);

    *trsf =
      sqrt(_PI_/(2.*l+5.))*Sp/(l+2.5) * l*(l+2.)/(2.*l+1.)/(2.*l+3.)
      -sqrt(_PI_/(2.*l+1.))*S/(l+0.5) * l/(2.*l+1.)*(l/(2.*l-1.)+(l+1.)/(2.*l+3.))
      +sqrt(_PI_/(2.*l-3.))*Sm/(l-1.5) * l*(l-1.)/(2.*l+1.)/(2.*l-1.);

  }

  else {

    class_stop(ptr->error_message,
               "Limber approximation has not been coded for the radial_type of index %d\n",
               radial_type);

  }

  return _SUCCESS_;

}

int transfer_limber_interpolate(
                                struct transfers * ptr,
                                double * tau0_minus_tau,
                                double * sources,
                                int tau_size,
                                double tau0_minus_tau_limber,
                                double * S
                                ){

  int index_tau;
  double dS,ddS;

  /** - find  bracketing indices.
      index_tau must be at least 1 (so that index_tau-1 is at least 0)
      and at most tau_size-2 (so that index_tau+1 is at most tau_size-1).
  */
  index_tau=1;
  while ((tau0_minus_tau[index_tau] > tau0_minus_tau_limber) && (index_tau<tau_size-2))
    index_tau++;

  /** - interpolate by fitting a polynomial of order two; get source
      and its first two derivatives. Note that we are not
      interpolating S, but the product S*(tau0-tau). Indeed this
      product is regular in tau=tau0, while S alone diverges for
      lensing. */

  /* the case where the last of the three point is the edge (tau0=tau) must be treated separately, see below */
  if (index_tau < tau_size-2) {

    class_call(array_interpolate_parabola(tau0_minus_tau[index_tau-1],
                                          tau0_minus_tau[index_tau],
                                          tau0_minus_tau[index_tau+1],
                                          tau0_minus_tau_limber,
                                          sources[index_tau-1]*tau0_minus_tau[index_tau-1],
                                          sources[index_tau]*tau0_minus_tau[index_tau],
                                          sources[index_tau+1]*tau0_minus_tau[index_tau+1],
                                          S,
                                          &dS,
                                          &ddS,
                                          ptr->error_message),
               ptr->error_message,
               ptr->error_message);

  }

  /* in this case, we have stored a zero for sources[index_k*tau_size+index_tau+1]. But we can use in very good approximation the fact that S*(tau0-tau) is constant near tau=tau0 and replace sources[index_k*tau_size+index_tau+1]*tau0_minus_tau[index_tau+1] by sources[index_k*tau_size+index_tau]*tau0_minus_tau[index_tau] */
  else {

    class_call(array_interpolate_parabola(tau0_minus_tau[index_tau-1],
                                          tau0_minus_tau[index_tau],
                                          tau0_minus_tau[index_tau+1],
                                          tau0_minus_tau_limber,
                                          sources[index_tau-1]*tau0_minus_tau[index_tau-1],
                                          sources[index_tau]*tau0_minus_tau[index_tau],
                                          sources[index_tau]*tau0_minus_tau[index_tau],
                                          S,
                                          &dS,
                                          &ddS,
                                          ptr->error_message),
               ptr->error_message,
               ptr->error_message);
  }

  return _SUCCESS_;

}

/**
 * This routine computes the transfer functions \f$ \Delta_l^{X} (k)
 * \f$) for each mode, initial condition, type, multipole l and
 * wavenumber k, by using the Limber approximation at ordet two, i.e
 * as a function of the source function and its first two derivatives
 * at a single value of tau
 *
 * @param ppt                   Input : pointer to perturbation structure
 * @param ptr                   Input : pointer to transfers structure
 * @param tau0                  Input : conformal time today
 * @param index_md            Input : index of mode
 * @param index_tt              Input : index of type
 * @param index_l               Input : index of multipole
 * @param index_k               Input : index of wavenumber
 * @param interpolated_sources  Input: array of interpolated sources
 * @param trsf                  Output: transfer function \f$ \Delta_l(k) \f$
 * @return the error status
 */

int transfer_limber2(
                     int tau_size,
                     struct transfers * ptr,
                     int index_md,
                     int index_k,
                     double l,
                     double k,
                     double * tau0_minus_tau,
                     double * sources,
                     radial_function_type radial_type,
                     double * trsf
                     ){

  /** Summary: */

  /** - define local variables */

  /* conformal time at which source must be computed */
  double tau0_minus_tau_limber;
  int index_tau;

  /* interpolated source and its derivatives */
  double S, dS, ddS;

  /** - get k, l and infer tau such that k(tau0-tau)=l+1/2;
      check that tau is in appropriate range */

  tau0_minus_tau_limber = (l+0.5)/k;  //TBC: to be updated to include curvature effects

  if ((tau0_minus_tau_limber > tau0_minus_tau[0]) ||
      (tau0_minus_tau_limber < tau0_minus_tau[tau_size-1])) {
    *trsf = 0.;
    return _SUCCESS_;
  }

  /** - find  bracketing indices */
  index_tau=0;
  while ((tau0_minus_tau[index_tau] > tau0_minus_tau_limber) && (index_tau<tau_size-2))
    index_tau++;

  /** - interpolate by fitting a polynomial of order two; get source
      and its first two derivatives */
  class_call(array_interpolate_parabola(tau0_minus_tau[index_tau-1],
                                        tau0_minus_tau[index_tau],
                                        tau0_minus_tau[index_tau+1],
                                        tau0_minus_tau_limber,
                                        sources[index_tau-1],
                                        sources[index_tau],
                                        sources[index_tau+1],
                                        &S,
                                        &dS,
                                        &ddS,
                                        ptr->error_message),
             ptr->error_message,
             ptr->error_message);


  /** - get transfer from 2nd order Limber approx (infered from 0809.5112 [astro-ph]) */

  *trsf = sqrt(_PI_/(2.*l+1.))/k*((1.-3./2./(2.*l+1.)/(2.*l+1.))*S+dS/k/(2.*l+1.)-0.5*ddS/k/k);

  return _SUCCESS_;

}

int transfer_can_be_neglected(
                              struct precision * ppr,
                              struct perturbs * ppt,
                              struct transfers * ptr,
                              int index_md,
                              int index_ic,
                              int index_tt,
                              double ra_rec,
                              double k,
                              double l,
                              short * neglect) {

  *neglect = _FALSE_;

  if (_scalars_) {

    if ((ppt->has_cl_cmb_temperature == _TRUE_) && (index_tt == ptr->index_tt_t0) && (l < (k-ppr->transfer_neglect_delta_k_S_t0)*ra_rec)) *neglect = _TRUE_;

    else if ((ppt->has_cl_cmb_temperature == _TRUE_) && (index_tt == ptr->index_tt_t1) && (l < (k-ppr->transfer_neglect_delta_k_S_t1)*ra_rec)) *neglect = _TRUE_;

    else if ((ppt->has_cl_cmb_temperature == _TRUE_) && (index_tt == ptr->index_tt_t2) && (l < (k-ppr->transfer_neglect_delta_k_S_t2)*ra_rec)) *neglect = _TRUE_;

    else if ((ppt->has_cl_cmb_polarization == _TRUE_) && (index_tt == ptr->index_tt_e) && (l < (k-ppr->transfer_neglect_delta_k_S_e)*ra_rec)) *neglect = _TRUE_;

  }

  else if (_vectors_) {

    if ((ppt->has_cl_cmb_temperature == _TRUE_) && (index_tt == ptr->index_tt_t1) && (l < (k-ppr->transfer_neglect_delta_k_V_t1)*ra_rec)) *neglect = _TRUE_;

    else if ((ppt->has_cl_cmb_temperature == _TRUE_) && (index_tt == ptr->index_tt_t2) && (l < (k-ppr->transfer_neglect_delta_k_V_t2)*ra_rec)) *neglect = _TRUE_;

    else if ((ppt->has_cl_cmb_polarization == _TRUE_) && (index_tt == ptr->index_tt_e) && (l < (k-ppr->transfer_neglect_delta_k_V_e)*ra_rec)) *neglect = _TRUE_;

    else if ((ppt->has_cl_cmb_polarization == _TRUE_) && (index_tt == ptr->index_tt_b) && (l < (k-ppr->transfer_neglect_delta_k_V_b)*ra_rec)) *neglect = _TRUE_;

  }

  else if (_tensors_) {

    if ((ppt->has_cl_cmb_temperature == _TRUE_) && (index_tt == ptr->index_tt_t2) && (l < (k-ppr->transfer_neglect_delta_k_T_t2)*ra_rec)) *neglect = _TRUE_;

    else if ((ppt->has_cl_cmb_polarization == _TRUE_) && (index_tt == ptr->index_tt_e) && (l < (k-ppr->transfer_neglect_delta_k_T_e)*ra_rec)) *neglect = _TRUE_;

    else if ((ppt->has_cl_cmb_polarization == _TRUE_) && (index_tt == ptr->index_tt_b) && (l < (k-ppr->transfer_neglect_delta_k_T_b)*ra_rec)) *neglect = _TRUE_;

  }

  return _SUCCESS_;

}

int transfer_late_source_can_be_neglected(
                                          struct precision * ppr,
                                          struct perturbs * ppt,
                                          struct transfers * ptr,
                                          int index_md,
                                          int index_tt,
                                          double l,
                                          short * neglect) {

  *neglect = _FALSE_;

  if (l > ppr->transfer_neglect_late_source*ptr->angular_rescaling) {

    /* sources at late times canb be neglected for CMB, excepted when
       there is a LISW: this means for tt_t1, t2, e */

    if (_scalars_) {
      if (ppt->has_cl_cmb_temperature == _TRUE_) {
        if ((index_tt == ptr->index_tt_t1) ||
            (index_tt == ptr->index_tt_t2))
          *neglect = _TRUE_;
      }
      if (ppt->has_cl_cmb_polarization == _TRUE_) {
        if (index_tt == ptr->index_tt_e)
          *neglect = _TRUE_;
      }
    }
    else if (_vectors_) {
      if (ppt->has_cl_cmb_temperature == _TRUE_) {
        if ((index_tt == ptr->index_tt_t1) ||
            (index_tt == ptr->index_tt_t2))
          *neglect = _TRUE_;
      }
      if (ppt->has_cl_cmb_polarization == _TRUE_) {
        if ((index_tt == ptr->index_tt_e) ||
            (index_tt == ptr->index_tt_b))
          *neglect = _TRUE_;
      }
    }
    else if (_tensors_) {
      if (ppt->has_cl_cmb_polarization == _TRUE_) {
        if ((index_tt == ptr->index_tt_e) ||
            (index_tt == ptr->index_tt_b))
          *neglect = _TRUE_;
      }
    }
  }

  return _SUCCESS_;

}

int transfer_radial_function(
                             struct transfer_workspace * ptw,
                             struct perturbs * ppt,
                             struct transfers * ptr,
                             double k,
                             int index_q,
                             int index_l,
                             int x_size,
                             double * radial_function,
                             radial_function_type radial_type
                             ){

  HyperInterpStruct * pHIS;
  double *chi = ptw->chi;
  double *cscKgen = ptw->cscKgen;
  double *cotKgen = ptw->cotKgen;
  int j;
  double *Phi, *dPhi, *d2Phi, *chireverse;
  double K=0.,k2=1.0;
  double sqrt_absK_over_k;
  double absK_over_k2;
  double nu=0., chi_tp=0.;
  double factor, s0, s2, ssqrt3, si, ssqrt2, ssqrt2i;
  double l = (double)ptr->l[index_l];
  double rescale_argument;
  double rescale_amplitude;
  double * rescale_function;
  int (*interpolate_Phi)();
  int (*interpolate_dPhi)();
  int (*interpolate_Phid2Phi)();
  int (*interpolate_PhidPhi)();
  int (*interpolate_PhidPhid2Phi)();
  enum Hermite_Interpolation_Order HIorder;

  K = ptw->K;
  k2 = k*k;

  if (ptw->sgnK==0){
    /* This is the choice consistent with chi=k*(tau0-tau) and nu=1 */
    sqrt_absK_over_k = 1.0;
  }
  else {
    K=ptw->K;
    sqrt_absK_over_k = sqrt(ptw->sgnK*K)/k;
  }
  absK_over_k2 =sqrt_absK_over_k*sqrt_absK_over_k;

  class_alloc(Phi,sizeof(double)*x_size,ptr->error_message);
  class_alloc(dPhi,sizeof(double)*x_size,ptr->error_message);
  class_alloc(d2Phi,sizeof(double)*x_size,ptr->error_message);
  class_alloc(chireverse,sizeof(double)*x_size,ptr->error_message);
  class_alloc(rescale_function,sizeof(double)*x_size,ptr->error_message);

  if (ptw->sgnK == 0) {
    pHIS = ptw->pBIS;
    rescale_argument = 1.;
    rescale_amplitude = 1.;
    HIorder = HERMITE4;
  }
  else if (index_q < ptr->index_q_flat_approximation) {
    pHIS = &(ptw->HIS);
    rescale_argument = 1.;
    rescale_amplitude = 1.;
    HIorder = HERMITE6;
  }
  else {
    pHIS = ptw->pBIS;
    if (ptw->sgnK == 1){
      nu = ptr->q[index_q]/sqrt(K);
      chi_tp = asin(sqrt(ptr->l[index_l]*(ptr->l[index_l]+1.))/nu);
    }
    else{
      nu = ptr->q[index_q]/sqrt(-K);
      chi_tp = asinh(sqrt(ptr->l[index_l]*(ptr->l[index_l]+1.))/nu);
    }
    rescale_argument = sqrt(ptr->l[index_l]*(ptr->l[index_l]+1.))/chi_tp;
    rescale_amplitude = pow(1.-K*ptr->l[index_l]*(ptr->l[index_l]+1.)/ptr->q[index_q]/ptr->q[index_q],-1./12.);
    HIorder = HERMITE4;
  }

  switch (HIorder){
  case HERMITE3:
    interpolate_Phi = hyperspherical_Hermite3_interpolation_vector_Phi;
    interpolate_dPhi = hyperspherical_Hermite3_interpolation_vector_dPhi;
    interpolate_PhidPhi = hyperspherical_Hermite3_interpolation_vector_PhidPhi;
    interpolate_Phid2Phi = hyperspherical_Hermite3_interpolation_vector_Phid2Phi;
    interpolate_PhidPhid2Phi = hyperspherical_Hermite3_interpolation_vector_PhidPhid2Phi;
    break;
  case HERMITE4:
    interpolate_Phi = hyperspherical_Hermite4_interpolation_vector_Phi;
    interpolate_dPhi = hyperspherical_Hermite4_interpolation_vector_dPhi;
    interpolate_PhidPhi = hyperspherical_Hermite4_interpolation_vector_PhidPhi;
    interpolate_Phid2Phi = hyperspherical_Hermite4_interpolation_vector_Phid2Phi;
    interpolate_PhidPhid2Phi = hyperspherical_Hermite4_interpolation_vector_PhidPhid2Phi;
    break;
  case HERMITE6:
    interpolate_Phi = hyperspherical_Hermite6_interpolation_vector_Phi;
    interpolate_dPhi = hyperspherical_Hermite6_interpolation_vector_dPhi;
    interpolate_PhidPhi = hyperspherical_Hermite6_interpolation_vector_PhidPhi;
    interpolate_Phid2Phi = hyperspherical_Hermite6_interpolation_vector_Phid2Phi;
    interpolate_PhidPhid2Phi = hyperspherical_Hermite6_interpolation_vector_PhidPhid2Phi;
    break;
  }

  //Reverse chi
  for (j=0; j<x_size; j++) {
    chireverse[j] = chi[x_size-1-j]*rescale_argument;
    if (rescale_amplitude == 1.) {
      rescale_function[j] = 1.;
    }
    else {
      if (ptw->sgnK == 1) {
        rescale_function[j] =
          MIN(
              rescale_amplitude
              * (1
                 + 0.34 * atan(ptr->l[index_l]/nu) * (chireverse[j]/rescale_argument-chi_tp)
                 + 2.00 * pow(atan(ptr->l[index_l]/nu) * (chireverse[j]/rescale_argument-chi_tp),2)),
              chireverse[j]/rescale_argument/sin(chireverse[j]/rescale_argument)
              );
      }
      else {
        rescale_function[j] =
          MAX(
              rescale_amplitude
              * (1
                 - 0.38 * atan(ptr->l[index_l]/nu) * (chireverse[j]/rescale_argument-chi_tp)
                 + 0.40 * pow(atan(ptr->l[index_l]/nu) * (chireverse[j]/rescale_argument-chi_tp),2)),
              chireverse[j]/rescale_argument/sinh(chireverse[j]/rescale_argument)
              );
      }
    }
  }

  /*
    class_test(pHIS->x[0] > chireverse[0],
    ptr->error_message,
    "Bessels need to be interpolated at %e, outside the range in which they have been computed (>%e). Decrease their x_min.",
    chireverse[0],
    pHIS->x[0]);
  */

  class_test((pHIS->x[pHIS->x_size-1] < chireverse[x_size-1]) && (ptw->sgnK != 1),
             ptr->error_message,
             "Bessels need to be interpolated at %e, outside the range in which they have been computed (<%e). Increase their x_max.",
             chireverse[x_size-1],
             pHIS->x[pHIS->x_size-1]
             );

  switch (radial_type){
  case SCALAR_TEMPERATURE_0:
    class_call(interpolate_Phi(pHIS, x_size, index_l, chireverse, Phi, ptr->error_message),
               ptr->error_message, ptr->error_message);
    //hyperspherical_Hermite_interpolation_vector(pHIS, x_size, index_l, chireverse, Phi, NULL, NULL);
    for (j=0; j<x_size; j++)
      radial_function[x_size-1-j] = Phi[j]*rescale_function[j];
    break;
  case SCALAR_TEMPERATURE_1:
    class_call(interpolate_dPhi(pHIS, x_size, index_l, chireverse, dPhi, ptr->error_message),
               ptr->error_message, ptr->error_message);
    //hyperspherical_Hermite_interpolation_vector(pHIS, x_size, index_l, chireverse, NULL, dPhi, NULL);
    for (j=0; j<x_size; j++)
      radial_function[x_size-1-j] = sqrt_absK_over_k*dPhi[j]*rescale_argument*rescale_function[j];
    break;
  case SCALAR_TEMPERATURE_2:
    class_call(interpolate_Phid2Phi(pHIS, x_size, index_l, chireverse, Phi, d2Phi, ptr->error_message),
               ptr->error_message, ptr->error_message);
    //hyperspherical_Hermite_interpolation_vector(pHIS, x_size, index_l, chireverse, Phi, NULL, d2Phi);
    s2 = sqrt(1.0-3.0*K/k2);
    factor = 1.0/(2.0*s2);
    for (j=0; j<x_size; j++)
      radial_function[x_size-1-j] = factor*(3*absK_over_k2*d2Phi[j]*rescale_argument*rescale_argument+Phi[j])*rescale_function[j];
    break;
  case SCALAR_POLARISATION_E:
    class_call(interpolate_Phi(pHIS, x_size, index_l, chireverse, Phi, ptr->error_message),
               ptr->error_message, ptr->error_message);
    //hyperspherical_Hermite_interpolation_vector(pHIS, x_size, index_l, chireverse, Phi, NULL, NULL);
    s2 = sqrt(1.0-3.0*K/k2);
    factor = sqrt(3.0/8.0*(l+2.0)*(l+1.0)*l*(l-1.0))/s2;
    for (j=0; j<x_size; j++)
      radial_function[x_size-1-j] = factor*cscKgen[x_size-1-j]*cscKgen[x_size-1-j]*Phi[j]*rescale_function[j];
    break;
  case VECTOR_TEMPERATURE_1:
    class_call(interpolate_Phi(pHIS, x_size, index_l, chireverse, Phi, ptr->error_message),
               ptr->error_message, ptr->error_message);
    //hyperspherical_Hermite_interpolation_vector(pHIS, x_size, index_l, chireverse, Phi, NULL, NULL);
    s0 = sqrt(1.0+K/k2);
    factor = sqrt(0.5*l*(l+1))/s0;
    for (j=0; j<x_size; j++)
      radial_function[x_size-1-j] = factor*cscKgen[x_size-1-j]*Phi[j]*rescale_function[j];
    break;
  case VECTOR_TEMPERATURE_2:
    class_call(interpolate_PhidPhi(pHIS, x_size, index_l, chireverse, Phi, dPhi, ptr->error_message),
               ptr->error_message, ptr->error_message);
    //hyperspherical_Hermite_interpolation_vector(pHIS, x_size, index_l, chireverse, Phi, dPhi, NULL);
    s0 = sqrt(1.0+K/k2);
    ssqrt3 = sqrt(1.0-2.0*K/k2);
    factor = sqrt(1.5*l*(l+1))/s0/ssqrt3;
    for (j=0; j<x_size; j++)
      radial_function[x_size-1-j] = factor*cscKgen[x_size-1-j]*(sqrt_absK_over_k*dPhi[j]*rescale_argument-cotKgen[j]*Phi[j])*rescale_function[j];
    break;
  case VECTOR_POLARISATION_E:
    class_call(interpolate_PhidPhi(pHIS, x_size, index_l, chireverse, Phi, dPhi, ptr->error_message),
               ptr->error_message, ptr->error_message);
    //    hyperspherical_Hermite_interpolation_vector(pHIS, x_size, index_l, chireverse, Phi, dPhi, NULL);
    s0 = sqrt(1.0+K/k2);
    ssqrt3 = sqrt(1.0-2.0*K/k2);
    factor = 0.5*sqrt((l-1.0)*(l+2.0))/s0/ssqrt3;
    for (j=0; j<x_size; j++)
      radial_function[x_size-1-j] = factor*cscKgen[x_size-1-j]*(cotKgen[j]*Phi[j]+sqrt_absK_over_k*dPhi[j]*rescale_argument)*rescale_function[j];
    break;
  case VECTOR_POLARISATION_B:
    class_call(interpolate_Phi(pHIS, x_size, index_l, chireverse, Phi, ptr->error_message),
               ptr->error_message, ptr->error_message);
    //hyperspherical_Hermite_interpolation_vector(pHIS, x_size, index_l, chireverse, Phi, NULL, NULL);
    s0 = sqrt(1.0+K/k2);
    ssqrt3 = sqrt(1.0-2.0*K/k2);
    si = sqrt(1.0+2.0*K/k2);
    factor = 0.5*sqrt((l-1.0)*(l+2.0))*si/s0/ssqrt3;
    for (j=0; j<x_size; j++)
      radial_function[x_size-1-j] = factor*cscKgen[x_size-1-j]*Phi[j]*rescale_function[j];
    break;
  case TENSOR_TEMPERATURE_2:
    class_call(interpolate_Phi(pHIS, x_size, index_l, chireverse, Phi, ptr->error_message),
               ptr->error_message, ptr->error_message);
    //hyperspherical_Hermite_interpolation_vector(pHIS, x_size, index_l, chireverse, Phi, NULL, NULL);
    ssqrt2 = sqrt(1.0-1.0*K/k2);
    si = sqrt(1.0+2.0*K/k2);
    factor = sqrt(3.0/8.0*(l+2.0)*(l+1.0)*l*(l-1.0))/si/ssqrt2;
    for (j=0; j<x_size; j++)
      radial_function[x_size-1-j] = factor*cscKgen[x_size-1-j]*cscKgen[x_size-1-j]*Phi[j]*rescale_function[j];
    break;
  case TENSOR_POLARISATION_E:
    class_call(interpolate_PhidPhid2Phi(pHIS, x_size, index_l, chireverse, Phi, dPhi, d2Phi, ptr->error_message),
               ptr->error_message, ptr->error_message);
    //hyperspherical_Hermite_interpolation_vector(pHIS, x_size, index_l, chireverse, Phi, NULL, NULL);
    ssqrt2 = sqrt(1.0-1.0*K/k2);
    si = sqrt(1.0+2.0*K/k2);
    factor = 0.25/si/ssqrt2;
    for (j=0; j<x_size; j++)
      radial_function[x_size-1-j] = factor*(absK_over_k2*d2Phi[j]*rescale_argument*rescale_argument
                                            +4.0*cotKgen[x_size-1-j]*sqrt_absK_over_k*dPhi[j]*rescale_argument
                                            -(1.0+4*K/k2-2.0*cotKgen[x_size-1-j]*cotKgen[x_size-1-j])*Phi[j])*rescale_function[j];
    break;
  case TENSOR_POLARISATION_B:
    class_call(interpolate_PhidPhi(pHIS, x_size, index_l, chireverse, Phi, dPhi, ptr->error_message),
               ptr->error_message, ptr->error_message);
    //hyperspherical_Hermite_interpolation_vector(pHIS, x_size, index_l, chireverse, Phi, dPhi, NULL);
    ssqrt2i = sqrt(1.0+3.0*K/k2);
    ssqrt2 = sqrt(1.0-1.0*K/k2);
    si = sqrt(1.0+2.0*K/k2);
    factor = 0.5*ssqrt2i/ssqrt2/si;
    for (j=0; j<x_size; j++)
      radial_function[x_size-1-j] = factor*(sqrt_absK_over_k*dPhi[j]*rescale_argument+2.0*cotKgen[x_size-1-j]*Phi[j])*rescale_function[j];
    break;
  case NC_RSD:
    class_call(interpolate_Phid2Phi(pHIS, x_size, index_l, chireverse, Phi, d2Phi, ptr->error_message),
               ptr->error_message, ptr->error_message);
    //hyperspherical_Hermite_interpolation_vector(pHIS, x_size, index_l, chireverse, Phi, NULL, d2Phi);
    //s2 = sqrt(1.0-3.0*K/k2);
    factor = 1.0;
    for (j=0; j<x_size; j++)
      radial_function[x_size-1-j] = factor*absK_over_k2*d2Phi[j]*rescale_argument*rescale_argument*rescale_function[j];
      // Note: in previous line there was a missing factor absK_over_k2 until version 2.4.3. Credits Francesco Montanari.
    break;
  }

  free(Phi);
  free(dPhi);
  free(d2Phi);
  free(chireverse);
  free(rescale_function);

  return _SUCCESS_;
}

int transfer_select_radial_function(
                                    struct perturbs * ppt,
                                    struct transfers * ptr,
                                    int index_md,
                                    int index_tt,
                                    radial_function_type * radial_type
                                    ) {

  /* generic case leading to generic bessel function (it applies also to all nonCMB types: lcmb, density, lensing) */
  *radial_type = SCALAR_TEMPERATURE_0;

  /* other specific cases */
  if (_scalars_) {

    if (ppt->has_cl_cmb_temperature == _TRUE_) {

      if (index_tt == ptr->index_tt_t0) {
        *radial_type = SCALAR_TEMPERATURE_0;
      }
      if (index_tt == ptr->index_tt_t1) {
        *radial_type = SCALAR_TEMPERATURE_1;
      }
      if (index_tt == ptr->index_tt_t2) {
        *radial_type = SCALAR_TEMPERATURE_2;
      }

    }

    if (ppt->has_cl_cmb_polarization == _TRUE_) {

      if (index_tt == ptr->index_tt_e) {
        *radial_type = SCALAR_POLARISATION_E;
      }

    }

    if (_index_tt_in_range_(ptr->index_tt_d1,      ptr->selection_num_nc, ppt->has_nc_rsd3))
      *radial_type = SCALAR_TEMPERATURE_1;

    if (_index_tt_in_range_(ptr->index_tt_rsd,     ptr->selection_num_nc, ppt->has_nc_rsd1))
      *radial_type = NC_RSD;

    if (_index_tt_in_range_(ptr->index_tt_nc_g5,   ptr->selection_num_nc, ppt->has_nc_gr5))
      *radial_type = SCALAR_TEMPERATURE_1;
  }

  if (_vectors_) {

    if (ppt->has_cl_cmb_temperature == _TRUE_) {

      if (index_tt == ptr->index_tt_t1) {
        *radial_type = VECTOR_TEMPERATURE_1;
      }
      if (index_tt == ptr->index_tt_t2) {
        *radial_type = VECTOR_TEMPERATURE_2;
      }
    }

    if (ppt->has_cl_cmb_polarization == _TRUE_) {

      if (index_tt == ptr->index_tt_e) {
        *radial_type = VECTOR_POLARISATION_E;
      }
      if (index_tt == ptr->index_tt_b) {
        *radial_type = VECTOR_POLARISATION_B;
      }

    }
  }

  if (_tensors_) {

    if (ppt->has_cl_cmb_temperature == _TRUE_) {

      if (index_tt == ptr->index_tt_t2) {
        *radial_type = TENSOR_TEMPERATURE_2;
      }
    }

    if (ppt->has_cl_cmb_polarization == _TRUE_) {

      if (index_tt == ptr->index_tt_e) {
        *radial_type = TENSOR_POLARISATION_E;
      }
      if (index_tt == ptr->index_tt_b) {
        *radial_type = TENSOR_POLARISATION_B;
      }

    }
  }

  return _SUCCESS_;

}

/* for reading global selection function (ie the one multiplying the selection function of each bin) */

int transfer_global_selection_read(
				   struct perturbs * ppt,
                                   struct transfers * ptr
                                   ) {
  int ii;
  FILE * input_file;
  int row,status;
  double tmp1,tmp2;

  for(ii=0;ii<ptr->n_tracers_nc;ii++) {
    //Read selection function
    if((ppt->has_nc_density) || (ppt->has_nc_rsd1) || (ppt->has_nc_rsd2) || (ppt->has_nc_rsd3) ||
       (ppt->has_nc_gr1) || (ppt->has_nc_gr2) || (ppt->has_nc_gr3) || (ppt->has_nc_lens) ||
       (ppt->has_nc_gr4) || (ppt->has_nc_gr5)) {
      ptr->nz_size_nc[ii] = 0;
      
      input_file = fopen(ptr->nz_file_name_nc[ii],"r");
      class_test(input_file == NULL,
		 ptr->error_message,
		 "Could not open file %s!",ptr->nz_file_name_nc[ii]);
      
      /* Find size of table */
      for (row=0,status=2; status==2; row++){
	status = fscanf(input_file,"%lf %lf",&tmp1,&tmp2);
      }
      rewind(input_file);
      ptr->nz_size_nc[ii] = row-1;
      
      /* Allocate room for interpolation table */
      class_alloc(ptr->nz_z_nc[ii],sizeof(double)*ptr->nz_size_nc[ii],ptr->error_message);
      class_alloc(ptr->nz_nz_nc[ii],sizeof(double)*ptr->nz_size_nc[ii],ptr->error_message);
      class_alloc(ptr->nz_ddnz_nc[ii],sizeof(double)*ptr->nz_size_nc[ii],ptr->error_message);
      
      for (row=0; row<ptr->nz_size_nc[ii]; row++){
	status = fscanf(input_file,"%lf %lf",
			&ptr->nz_z_nc[ii][row],&ptr->nz_nz_nc[ii][row]);
      }
      fclose(input_file);
      
      /* Call spline interpolation: */
      class_call(array_spline_table_lines(ptr->nz_z_nc[ii],
					  ptr->nz_size_nc[ii],
					  ptr->nz_nz_nc[ii],
					  1,
					  ptr->nz_ddnz_nc[ii],
					  _SPLINE_EST_DERIV_,
					  ptr->error_message),
		 ptr->error_message,
		 ptr->error_message);
    }

    //Read clustering bias
    if(ppt->has_nc_density) {
#ifdef _CLASST_DEBUG
      printf("Reading bias from %s\n",ptr->bz_file_name[ii]);
#endif //_CLASST_DEBUG
      input_file = fopen(ptr->bz_file_name[ii],"r");
      class_test(input_file == NULL,
		 ptr->error_message,
		 "Could not open file %s!",ptr->bz_file_name[ii]);
      
      /* Find size of table */
      for (row=0,status=2; status==2; row++){
	status = fscanf(input_file,"%lf %lf",&tmp1,&tmp2);
      }
      rewind(input_file);
      ptr->bz_size[ii] = row-1;
      
      /* Allocate room for interpolation table */
      class_alloc(ptr->bz_z[ii],sizeof(double)*ptr->bz_size[ii],ptr->error_message);
      class_alloc(ptr->bz_bz[ii],sizeof(double)*ptr->bz_size[ii],ptr->error_message);
      class_alloc(ptr->bz_ddbz[ii],sizeof(double)*ptr->bz_size[ii],ptr->error_message);
      
      for (row=0; row<ptr->bz_size[ii]; row++){
	status = fscanf(input_file,"%lf %lf",
			&ptr->bz_z[ii][row],&ptr->bz_bz[ii][row]);
      }
      fclose(input_file);
      
      /* Call spline interpolation: */
      class_call(array_spline_table_lines(ptr->bz_z[ii],
					  ptr->bz_size[ii],
					  ptr->bz_bz[ii],
					  1,
					  ptr->bz_ddbz[ii],
					  _SPLINE_EST_DERIV_,
					  ptr->error_message),
		 ptr->error_message,
		 ptr->error_message);
    }

    //Read magnification bias
    if((ppt->has_nc_rsd3) || (ppt->has_nc_gr2) || (ppt->has_nc_lens) ||
       (ppt->has_nc_gr4) || (ppt->has_nc_gr5)) {
#ifdef _CLASST_DEBUG
      printf("Reading s-bias from %s\n",ptr->sz_file_name[ii]);
#endif //_CLASST_DEBUG
      input_file = fopen(ptr->sz_file_name[ii],"r");
      class_test(input_file == NULL,
		 ptr->error_message,
		 "Could not open file %s!",ptr->sz_file_name[ii]);
      
      /* Find size of table */
      for (row=0,status=2; status==2; row++){
	status = fscanf(input_file,"%lf %lf",&tmp1,&tmp2);
      }
      rewind(input_file);
      ptr->sz_size[ii] = row-1;
    
      /* Allocate room for interpolation table */
      class_alloc(ptr->sz_z[ii],sizeof(double)*ptr->sz_size[ii],ptr->error_message);
      class_alloc(ptr->sz_sz[ii],sizeof(double)*ptr->sz_size[ii],ptr->error_message);
      class_alloc(ptr->sz_ddsz[ii],sizeof(double)*ptr->sz_size[ii],ptr->error_message);
      
      for (row=0; row<ptr->sz_size[ii]; row++){
	status = fscanf(input_file,"%lf %lf",
			&ptr->sz_z[ii][row],&ptr->sz_sz[ii][row]);
      }
      fclose(input_file);
      
      /* Call spline interpolation: */
      class_call(array_spline_table_lines(ptr->sz_z[ii],
					  ptr->sz_size[ii],
					  ptr->sz_sz[ii],
					  1,
					  ptr->sz_ddsz[ii],
					  _SPLINE_EST_DERIV_,
					  ptr->error_message),
		 ptr->error_message,
		 ptr->error_message);
    }

    //Read evolution bias
    if((ppt->has_nc_rsd2) || (ppt->has_nc_rsd3) || 
       (ppt->has_nc_gr2) || (ppt->has_nc_gr5)) {
#ifdef _CLASST_DEBUG
      printf("Reading e-bias from %s\n",ptr->ez_file_name[ii]);
#endif //_CLASST_DEBUG
      input_file = fopen(ptr->ez_file_name[ii],"r");
      class_test(input_file == NULL,
		 ptr->error_message,
		 "Could not open file %s!",ptr->ez_file_name[ii]);
      
      /* Find size of table */
      for (row=0,status=2; status==2; row++){
	status = fscanf(input_file,"%lf %lf",&tmp1,&tmp2);
      }
      rewind(input_file);
      ptr->ez_size[ii] = row-1;
      
      /* Allocate room for interpolation table */
      class_alloc(ptr->ez_z[ii],sizeof(double)*ptr->ez_size[ii],ptr->error_message);
      class_alloc(ptr->ez_ez[ii],sizeof(double)*ptr->ez_size[ii],ptr->error_message);
      class_alloc(ptr->ez_ddez[ii],sizeof(double)*ptr->ez_size[ii],ptr->error_message);
      
      for (row=0; row<ptr->ez_size[ii]; row++){
	status = fscanf(input_file,"%lf %lf",
			&ptr->ez_z[ii][row],&ptr->ez_ez[ii][row]);
      }
      fclose(input_file);
      
      /* Call spline interpolation: */
      class_call(array_spline_table_lines(ptr->ez_z[ii],
					  ptr->ez_size[ii],
					  ptr->ez_ez[ii],
					  1,
					  ptr->ez_ddez[ii],
					  _SPLINE_EST_DERIV_,
					  ptr->error_message),
		 ptr->error_message,
		 ptr->error_message);
    }
  }

  for(ii=0;ii<ptr->n_tracers_wl;ii++) {
    if((ppt->has_lensing_shear) || (ppt->has_intrinsic_alignment)) {
      ptr->nz_size_wl[ii] = 0;
      
      input_file = fopen(ptr->nz_file_name_wl[ii],"r");
      class_test(input_file == NULL,
		 ptr->error_message,
		 "Could not open file %s!",ptr->nz_file_name_wl[ii]);
      
      /* Find size of table */
      for (row=0,status=2; status==2; row++){
	status = fscanf(input_file,"%lf %lf",&tmp1,&tmp2);
      }
      rewind(input_file);
      ptr->nz_size_wl[ii] = row-1;
  
      /* Allocate room for interpolation table */
      class_alloc(ptr->nz_z_wl[ii],sizeof(double)*ptr->nz_size_wl[ii],ptr->error_message);
      class_alloc(ptr->nz_nz_wl[ii],sizeof(double)*ptr->nz_size_wl[ii],ptr->error_message);
      class_alloc(ptr->nz_ddnz_wl[ii],sizeof(double)*ptr->nz_size_wl[ii],ptr->error_message);
  
      for (row=0; row<ptr->nz_size_wl[ii]; row++){
	status = fscanf(input_file,"%lf %lf",
			&ptr->nz_z_wl[ii][row],&ptr->nz_nz_wl[ii][row]);
      }
      fclose(input_file);
  
      /* Call spline interpolation: */
      class_call(array_spline_table_lines(ptr->nz_z_wl[ii],
					  ptr->nz_size_wl[ii],
					  ptr->nz_nz_wl[ii],
					  1,
					  ptr->nz_ddnz_wl[ii],
					  _SPLINE_EST_DERIV_,
					  ptr->error_message),
		 ptr->error_message,
		 ptr->error_message);
    }

    if(ppt->has_intrinsic_alignment) {
#ifdef _CLASST_DEBUG
      printf("Reading alignment bias from %s\n",ptr->aI_file_name[ii]);
#endif //_CLASST_DEBUG
      input_file = fopen(ptr->aI_file_name[ii],"r");
      class_test(input_file == NULL,
		 ptr->error_message,
		 "Could not open file %s!",ptr->aI_file_name[ii]);
      
      /* Find size of table */
      for (row=0,status=2; status==2; row++){
	status = fscanf(input_file,"%lf %lf",&tmp1,&tmp2);
      }
      rewind(input_file);
      ptr->aI_size[ii] = row-1;
      
      /* Allocate room for interpolation table */
      class_alloc(ptr->aI_z[ii],sizeof(double)*ptr->aI_size[ii],ptr->error_message);
      class_alloc(ptr->aI_aI[ii],sizeof(double)*ptr->aI_size[ii],ptr->error_message);
      class_alloc(ptr->aI_ddaI[ii],sizeof(double)*ptr->aI_size[ii],ptr->error_message);
      
      for (row=0; row<ptr->aI_size[ii]; row++){
	status = fscanf(input_file,"%lf %lf",
			&ptr->aI_z[ii][row],&ptr->aI_aI[ii][row]);
      }
      fclose(input_file);
      //      exit(1);
      
      /* Call spline interpolation: */
      class_call(array_spline_table_lines(ptr->aI_z[ii],
					  ptr->aI_size[ii],
					  ptr->aI_aI[ii],
					  1,
					  ptr->aI_ddaI[ii],
					  _SPLINE_EST_DERIV_,
					  ptr->error_message),
		 ptr->error_message,
		 ptr->error_message);

#ifdef _CLASST_DEBUG
      printf("Reading red fraction from %s\n",ptr->fred_file_name[ii]);
#endif //_CLASST_DEBUG
      input_file = fopen(ptr->fred_file_name[ii],"r");
      class_test(input_file == NULL,
		 ptr->error_message,
		 "Could not open file %s!",ptr->fred_file_name[ii]);
      
      /* Find size of table */
      for (row=0,status=2; status==2; row++){
	status = fscanf(input_file,"%lf %lf",&tmp1,&tmp2);
      }
      rewind(input_file);
      ptr->fred_size[ii] = row-1;
      
      /* Allocate room for interpolation table */
      class_alloc(ptr->fred_z[ii],sizeof(double)*ptr->fred_size[ii],ptr->error_message);
      class_alloc(ptr->fred_fred[ii],sizeof(double)*ptr->fred_size[ii],ptr->error_message);
      class_alloc(ptr->fred_ddfred[ii],sizeof(double)*ptr->fred_size[ii],ptr->error_message);
      
      for (row=0; row<ptr->fred_size[ii]; row++){
	status = fscanf(input_file,"%lf %lf",
			&ptr->fred_z[ii][row],&ptr->fred_fred[ii][row]);
      }
      fclose(input_file);
      
      /* Call spline interpolation: */
      class_call(array_spline_table_lines(ptr->fred_z[ii],
					  ptr->fred_size[ii],
					  ptr->fred_fred[ii],
					  1,
					  ptr->fred_ddfred[ii],
					  _SPLINE_EST_DERIV_,
					  ptr->error_message),
		 ptr->error_message,
		 ptr->error_message);
    }
  }

  return _SUCCESS_;
};

int transfer_workspace_init(
                            struct transfers * ptr,
                            struct precision * ppr,
                            struct transfer_workspace **ptw,
                            int perturb_tau_size,
                            int tau_size_max,
                            double K,
                            int sgnK,
                            double tau0_minus_tau_cut,
                            HyperInterpStruct * pBIS){

  class_calloc(*ptw,1,sizeof(struct transfer_workspace),ptr->error_message);

  (*ptw)->tau_size_max = tau_size_max;
  (*ptw)->l_size = ptr->l_size_max;
  (*ptw)->HIS_allocated=_FALSE_;
  (*ptw)->pBIS = pBIS;
  (*ptw)->K = K;
  (*ptw)->sgnK = sgnK;
  (*ptw)->tau0_minus_tau_cut = tau0_minus_tau_cut;
  (*ptw)->neglect_late_source = _FALSE_;

  class_alloc((*ptw)->interpolated_sources,perturb_tau_size*sizeof(double),ptr->error_message);
  class_alloc((*ptw)->sources,tau_size_max*sizeof(double),ptr->error_message);
  class_alloc((*ptw)->tau0_minus_tau,tau_size_max*sizeof(double),ptr->error_message);
  class_alloc((*ptw)->w_trapz,tau_size_max*sizeof(double),ptr->error_message);
  class_alloc((*ptw)->chi,tau_size_max*sizeof(double),ptr->error_message);
  class_alloc((*ptw)->cscKgen,tau_size_max*sizeof(double),ptr->error_message);
  class_alloc((*ptw)->cotKgen,tau_size_max*sizeof(double),ptr->error_message);

  return _SUCCESS_;
}

int transfer_workspace_free(
                            struct transfers * ptr,
                            struct transfer_workspace *ptw
                            ) {

  if (ptw->HIS_allocated==_TRUE_){
    //Free HIS structure:
    class_call(hyperspherical_HIS_free(&(ptw->HIS),ptr->error_message),
               ptr->error_message,
               ptr->error_message);
  }
  free(ptw->interpolated_sources);
  free(ptw->sources);
  free(ptw->tau0_minus_tau);
  free(ptw->w_trapz);
  free(ptw->chi);
  free(ptw->cscKgen);
  free(ptw->cotKgen);

  free(ptw);
  return _SUCCESS_;
}

int transfer_update_HIS(
                        struct precision * ppr,
                        struct transfers * ptr,
                        struct transfer_workspace * ptw,
                        int index_q,
                        double tau0
                        ) {

  double nu,new_nu;
  int int_nu;
  double xmin, xmax, sampling, phiminabs, xtol;
  double sqrt_absK;
  int l_size_max;
  int index_l_left,index_l_right;

  if (ptw->HIS_allocated == _TRUE_) {
    class_call(hyperspherical_HIS_free(&(ptw->HIS),ptr->error_message),
               ptr->error_message,
               ptr->error_message);
    ptw->HIS_allocated = _FALSE_;
  }

  if ((ptw->sgnK!=0) && (index_q < ptr->index_q_flat_approximation)) {

    xmin = ppr->hyper_x_min;

    sqrt_absK = sqrt(ptw->sgnK*ptw->K);

    xmax = sqrt_absK*tau0;
    nu = ptr->q[index_q]/sqrt_absK;

    if (ptw->sgnK == 1) {
      xmax = MIN(xmax,_PI_/2.0-ppr->hyper_x_min); //We only need solution on [0;pi/2]

      int_nu = (int)(nu+0.2);
      new_nu = (double)int_nu;
      class_test(nu-new_nu > 1.e-6,
                 ptr->error_message,
                 "problem in q list definition in closed case for index_q=%d, nu=%e, nu-int(nu)=%e",index_q,nu,nu-new_nu);
      nu = new_nu;

    }

    if (nu > ppr->hyper_nu_sampling_step)
      sampling = ppr->hyper_sampling_curved_high_nu;
    else
      sampling = ppr->hyper_sampling_curved_low_nu;

    /* find the highest value of l such that x_nonzero < xmax = sqrt(|K|) tau0. That will be l_max. */
    l_size_max = ptr->l_size_max;
    if (ptw->sgnK == 1)
      while ((double)ptr->l[l_size_max-1] >= nu)
        l_size_max--;

    if (ptw->sgnK == -1){
      xtol = ppr->hyper_x_tol;
      phiminabs = ppr->hyper_phi_min_abs;

      /** First try to find lmax using fast approximation: */
      index_l_left=0;
      index_l_right=l_size_max-1;
      class_call(transfer_get_lmax(hyperspherical_get_xmin_from_approx,
                                   ptw->sgnK,
                                   nu,
                                   ptr->l,
                                   l_size_max,
                                   phiminabs,
                                   xmax,
                                   xtol,
                                   &index_l_left,
                                   &index_l_right,
                                   ptr->error_message),
                 ptr->error_message,
                 ptr->error_message);

      /** Now use WKB approximation to eventually modify borders: */
      class_call(transfer_get_lmax(hyperspherical_get_xmin_from_Airy,
                                   ptw->sgnK,
                                   nu,
                                   ptr->l,
                                   l_size_max,
                                   phiminabs,
                                   xmax,
                                   xtol,
                                   &index_l_left,
                                   &index_l_right,
                                   ptr->error_message),
                 ptr->error_message,
                 ptr->error_message);
      l_size_max = index_l_right+1;
    }

    class_test(nu <= 0.,
               ptr->error_message,
               "nu=%e when index_q=%d, q=%e, K=%e, sqrt(|K|)=%e; instead nu should always be strictly positive",
               nu,index_q,ptr->q[index_q],ptw->K,sqrt_absK);

    class_call(hyperspherical_HIS_create(ptw->sgnK,
                                         nu,
                                         l_size_max,
                                         ptr->l,
                                         xmin,
                                         xmax,
                                         sampling,
                                         ptr->l[l_size_max-1]+1,
                                         ppr->hyper_phi_min_abs,
                                         &(ptw->HIS),
                                         ptr->error_message),
               ptr->error_message,
               ptr->error_message);

    ptw->HIS_allocated = _TRUE_;

  }

  return _SUCCESS_;
}

int transfer_get_lmax(int (*get_xmin_generic)(int sgnK,
                                              int l,
                                              double nu,
                                              double xtol,
                                              double phiminabs,
                                              double *x_nonzero,
                                              int *fevals),
                      int sgnK,
                      double nu,
                      int *lvec,
                      int lsize,
                      double phiminabs,
                      double xmax,
                      double xtol,
                      int *index_l_left,
                      int *index_l_right,
                      ErrorMsg error_message){
  double x_nonzero;
  int fevals=0, index_l_mid;
  int multiplier;
  int right_boundary_checked = _FALSE_;
  int hil=0,hir=0,bini=0;
  class_call(get_xmin_generic(sgnK,
                              lvec[0],
                              nu,
                              xtol,
                              phiminabs,
                              &x_nonzero,
                              &fevals),
             error_message,
             error_message);
  if (x_nonzero >= xmax){
    //x at left boundary is already larger than xmax.
    *index_l_right = MAX(lsize-1,1);
    return _SUCCESS_;
  }
  class_call(get_xmin_generic(sgnK,
                              lvec[lsize-1],
                              nu,
                              xtol,
                              phiminabs,
                              &x_nonzero,
                              &fevals),
             error_message,
             error_message);

  if (x_nonzero < xmax){
    //All Bessels are relevant
    *index_l_left = MAX(0,(lsize-2));
    return _SUCCESS_;
  }
  /** Hunt for left boundary: */
  for (multiplier=1; ;multiplier *= 5){
    hil++;
    class_call(get_xmin_generic(sgnK,
                                lvec[*index_l_left],
                                nu,
                                xtol,
                                phiminabs,
                                &x_nonzero,
                                &fevals),
               error_message,
               error_message);
    if (x_nonzero <= xmax){
      //Boundary found
      break;
    }
    else{
      //We can use current index_l_left as index_l_right:
      *index_l_right = *index_l_left;
      right_boundary_checked = _TRUE_;
    }
    //Update index_l_left:
    *index_l_left = (*index_l_left)-multiplier;
    if (*index_l_left<=0){
      *index_l_left = 0;
      break;
    }
  }
  /** If not found, hunt for right boundary: */
  if (right_boundary_checked == _FALSE_){
    for (multiplier=1; ;multiplier *= 5){
      hir++;
      class_call(get_xmin_generic(sgnK,
                                  lvec[*index_l_right],
                                  nu,
                                  xtol,
                                  phiminabs,
                                  &x_nonzero,
                                  &fevals),
                 error_message,
                 error_message);
      if (x_nonzero >= xmax){
        //Boundary found
        break;
      }
      else{
        //We can use current index_l_right as index_l_left:
        *index_l_left = *index_l_right;
      }
      //Update index_l_right:
      *index_l_right = (*index_l_right)+multiplier;
      if (*index_l_right>=(lsize-1)){
        *index_l_right = lsize-1;
        break;
      }
    }
  }
  //  int fevalshunt=fevals;
  fevals=0;
  //Do binary search
  while (((*index_l_right) - (*index_l_left)) > 1) {
    bini++;
    index_l_mid= (int)(0.5*((*index_l_right)+(*index_l_left)));
    class_call(get_xmin_generic(sgnK,
                                lvec[index_l_mid],
                                nu,
                                xtol,
                                phiminabs,
                                &x_nonzero,
                                &fevals),
               error_message,
               error_message);
    if (x_nonzero < xmax)
      *index_l_left=index_l_mid;
    else
      *index_l_right=index_l_mid;
  }
  return _SUCCESS_;
}
