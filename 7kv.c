#include "udf.h"
#include "mem.h"
#include "sg_udms.h"
#define q 1.0
#define eps 8.854E-12
#define current 3.818e-6
#define dens 1.205
#define b 0.00018
#define x_ionization -3.0e-3
#define factor (2.0 * M_PI)

DEFINE_SOURCE(ymom_source,c,t,dS,eqn)
	{
		real x[ND_ND];
		real source;
		real g_phi;
		g_phi=C_UDSI_G(c,t,0)[1];
		source = -q*C_UDSI(c,t,2)*g_phi;
		dS[eqn]=0;
		return source;
	}

DEFINE_SOURCE(xmom_source,c,t,dS,eqn)
	{
		real x[ND_ND];
		real source;
		real g_phi;
		g_phi=C_UDSI_G(c,t,0)[0];
		source = -q*C_UDSI(c,t,2)*g_phi;
		dS[eqn]=0;
		return source;
	}

DEFINE_SOURCE(phi_source,c,t,dS,eqn)
	{
		real x[ND_ND];
		real source;
		face_t f;
		source = dens * C_UDSI(c,t,2) / eps;
		dS[eqn]=0;
		return source;
	}


DEFINE_SOURCE(charge_source,c,t,dS,eqn)
	{
#if !RP_HOST
		real x[ND_ND];
		real source, volume1;
                real e;
                real e0;
		real e1;
		C_CENTROID(x, c, t);
		e = NV_MAG(C_UDSI_G(c,t,0));
		e0 = 3.24e6;
		e1 = 2.8e6;

		if (e<e0 && e>e1 && x[0] > x_ionization)
		{
			volume1 = C_UDMI(c, t, 0);
			source = dens * current / volume1;
			dS[eqn] = 0.0;
             return source;
		} 
		else {
			source = 0.0;
			dS[eqn] = 0.0;
			return source;
                }
#endif
	}


DEFINE_ADJUST(effective_volume,d)
{
	Thread *t;
	cell_t c;
    real e, e0, e1;
    real volume1;
	real x[ND_ND];
	e0 = 3.24e6;
	e1 = 2.8e6;
	volume1 = 1e-16;
	if (! Data_Valid_P())
        {
		return;
        }
	t = Lookup_Thread(d, 12);
	begin_c_loop_int(c,t)
	    {
	       C_CENTROID(x, c, t);
	       e = NV_MAG(C_UDSI_G(c,t,0));
	       if (e<e0 && e>e1 && x[0] > x_ionization)
		 {	
		   volume1 += C_VOLUME(c,t);
		 } 
	    }
	end_c_loop_int(c,t)

	volume1 = volume1 * factor;
	volume1 = PRF_GRSUM1(volume1);
	/*printf("Ionization volume: %g\n",volume1);*/
	begin_c_loop_int (c, t)
	{
		C_UDMI(c, t, 0) = volume1;
	}
	end_c_loop_int(c, t)
}

DEFINE_EXECUTE_AT_END(current_leak) /*flow is current leakage at outlet*/
{
#if !RP_HOST
	Domain *d;
	real flow = 0.;
	face_t f;
	Thread *t;
	d = Get_Domain(1);
	t = Lookup_Thread(d,16);
	begin_f_loop(f,t)
		if PRINCIPAL_FACE_P(f,t) /* test if the face is the principle face*/
	{
		flow+= F_FLUX(f,t)*F_UDSI(f,t,2);
	}
	end_f_loop(f,t)
	flow = flow * factor;
	flow = PRF_GRSUM1(flow);
	printf("Current leakage from outlet: %g\n",flow);
#endif
}

DEFINE_EXECUTE_AT_END(current_cathode) /*flow is current leakage at outlet*/
{
#if !RP_HOST
	Domain *d;
	real NV_VEC(A);
	real flow1 = 0.;
	cell_t c0 = -1;
	Thread *t0 = NULL;
	face_t f;
	Thread *t;
	d = Get_Domain(1);
	t = Lookup_Thread(d,14);
	begin_f_loop(f,t)
		if PRINCIPAL_FACE_P(f,t) /* test if the face is the principle face*/
	{
	F_AREA(A,f,t);	
	c0 = F_C0(f,t);
	t0 = F_C0_THREAD(f,t);
	flow1 += - b * F_UDSI(f,t,2) * NV_DOT(C_UDSI_G(c0,t0,0),A);
	}
	end_f_loop(f,t)
	flow1 = flow1 * factor;
	flow1 = PRF_GRSUM1(flow1);
	printf("Current collected from cathode: %g\n",flow1);
#endif
}

DEFINE_UDS_FLUX(uds1_flux,f,t,i)
	{
		real a, a0, a1;
		Domain *d;
		cell_t c0,c1=-1;
		Thread *t0,*t1 =NULL;
		real NV_VEC(psi_vec),NV_VEC(A),flux=0.0;
		c0=F_C0(f,t);
		t0=F_C0_THREAD(f,t);
		F_AREA(A,f,t);
		d = Get_Domain(1);
		if (BOUNDARY_FACE_THREAD_P(t))
		{
			psi_vec[0] = F_U(f,t) * dens;
			psi_vec[1] = F_V(f,t) * dens;
			a=- b * dens;
			NV_V_VS(psi_vec,=,psi_vec,+,C_UDSI_G(c0,t0,0),*,a);
			flux = NV_DOT(psi_vec,A);}
		else
		{
			c1=F_C1(f,t);
			t1=F_C1_THREAD(f,t);
			psi_vec[0] = (C_U(c0,t0) + C_U(c1,t1)) * dens;
			psi_vec[1] = (C_V(c0,t0) + C_V(c1,t1)) * dens;
			a0 = - b * C_R(c0,t0);
			a1 = - b * C_R(c1,t1);
			NV_V_VS(psi_vec,=,psi_vec,+,C_UDSI_G(c0,t0,0),*,a0);
			NV_V_VS(psi_vec,=,psi_vec,+,C_UDSI_G(c1,t1,0),*,a1);
			flux = NV_DOT(psi_vec,A) / 2.0;
		}
		if (t == Lookup_Thread(d,18))
			flux = 0.0;
		return flux;
	}

/*DEFINE_UDS_UNSTEADY(ION,c,t,i,apu,su)
	{
		real physical_dt, vol, rho, phi_old;
		physical_dt = RP_Get_Real("physical-time-step");
		vol = C_VOLUME(c,t);
		rho = C_R_M1(c,t);
		*apu = -rho*vol / physical_dt;
		phi_old = C_STORAGE_R(c,t,SV_UDSI_M1(2));
		*su = rho * vol * phi_old / physical_dt;
	}*/

DEFINE_ADJUST(electrical_field_x,d)
{
#if !RP_HOST
	Thread *t;
	face_t f;
	cell_t c;
	Thread *tf;
	int n;
	if (! Data_Valid_P())
        {
		return;
        }
	
	thread_loop_c(t,d)
	 {
	 if (FLUID_THREAD_P(t))
	  {
	  begin_c_loop(c,t)
	    {
			C_UDSI(c,t,1) = - C_UDSI_G(c,t,0)[0];
			c_face_loop(c,t,n)
			{
				f = C_FACE(c,t,n);
				tf = C_FACE_THREAD(c,t,n);
				if (BOUNDARY_FACE_THREAD_P(tf))
				{ if (THREAD_F_AXIS != THREAD_TYPE(tf))
				  F_UDSI(f,tf,1) = C_UDSI(c,t,1);
				}
			}
	    }
	  end_c_loop(c,t)
	 }
	}
#endif
}

DEFINE_ADJUST(electrical_field_y,d)
{
#if !RP_HOST
	Thread *t;
	face_t f;
	cell_t c;
	Thread *tf;
	int n;
	if (! Data_Valid_P())
        {
		return;
        }
	
	thread_loop_c(t,d)
	 {
	 if (FLUID_THREAD_P(t))
	  {
	  begin_c_loop(c,t)
	    {
			C_UDSI(c,t,4) = - C_UDSI_G(c,t,0)[1];
			c_face_loop(c,t,n)
			{
				f = C_FACE(c,t,n);
				tf = C_FACE_THREAD(c,t,n);
				if (BOUNDARY_FACE_THREAD_P(tf))
				{ if (THREAD_F_AXIS != THREAD_TYPE(tf))
				  F_UDSI(f,tf,1) = C_UDSI(c,t,1);
				}
			}
	    }
	  end_c_loop(c,t)
	 }
	}
#endif
}

DEFINE_ON_DEMAND(rename_UDvars)
{
	Set_User_Scalar_Name(0,"UDS0: Electric Potential");
	Set_User_Scalar_Name(1,"UDS1: Ex");
	Set_User_Scalar_Name(2,"UDS2: Charge Density");
	Set_User_Scalar_Name(3,"UDS3: Current Residual");
	Set_User_Scalar_Name(4,"UDS4: Ey");
}

DEFINE_ADJUST(cell_current_density,d)
{
#if !RP_HOST
	Thread *t;
	face_t f;
	cell_t c;
	Thread *tf;
	int n;
	real NV_VEC(vec),NV_VEC(A);
	if (! Data_Valid_P())
         {
         return;
         }
	
	thread_loop_c(t,d)
	{
	   if (FLUID_THREAD_P(t))
	      {
	      begin_c_loop(c,t)
	         {
	         C_UDSI(c,t,3) = 0.;
	         vec[0] = (C_U(c,t) - b * C_UDSI_G(c,t,0)[0]) * C_UDSI(c,t,2) + C_UDSI_DIFF(c,t,2) * C_UDSI_G(c,t,2)[0] / dens;
	         vec[1] = (C_V(c,t) - b * C_UDSI_G(c,t,0)[1]) * C_UDSI(c,t,2) + C_UDSI_DIFF(c,t,2) * C_UDSI_G(c,t,2)[1] / dens;
	         c_face_loop(c,t,n)
	            {
	            f = C_FACE(c,t,n);
	            tf = C_FACE_THREAD(c,t,n);
	            F_AREA(A,f,tf);
	            C_UDSI(c,t,3) += NV_DOT(vec,A);
	            if (BOUNDARY_FACE_THREAD_P(tf))
	               { if (THREAD_F_AXIS != THREAD_TYPE(tf))
	                  F_UDSI(f,tf,3) = C_UDSI(c,t,3) * factor;
		         }
	            }
	            C_UDSI(c,t,3) = C_UDSI(c,t,3) * factor;
	         }
	      end_c_loop(c,t)
	   }
	}
#endif
}

DEFINE_EXECUTE_AT_END(correct_charge) /*flow is current leakage at outlet*/
{
#if !RP_HOST
	Domain *d;
	cell_t c;
	Thread *t;
	real x[ND_ND];
	d = Get_Domain(1);
	if (! Data_Valid_P())
        {
		return;
        }
	thread_loop_c(t,d)
	 {
	 if (FLUID_THREAD_P(t))
	  {
	  begin_c_loop_int(c,t)
	    {
			C_CENTROID(x, c, t);
			if (C_UDSI(c, t, 2) < 0. || x[0] < x_ionization - 2.0e-3)
				C_UDSI(c, t, 2) = 0.;
	    }
	  end_c_loop_int(c,t)
	 }
	}

#endif
}