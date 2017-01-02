import numpy as np
import os as os
import sys as sys
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import copy

FS=16

class ParamFisher:
    """ Fisher matrix parameter """
    val=0.0
    dval=0.0
    name="str"
    label="$x$"
    isfree=False
    do_plot=True
    prior=-1

    def __init__(self,val,dval,name,label,isfree,do_plot,prior):
        self.val=val
        self.dval=dval
        self.name=name
        self.label=label
        self.isfree=isfree
        self.do_plot=do_plot
        self.prior=prior
        
def find_param(param_list,name):
    index=0
    for par in param_list:
        if par.name==name :
            return index
        index+=1
    print "No parameter "+name
    sys.exit()

def plot_fisher_single(params,name,fishermat,ax,fc,lw,ls,lc,fact_axis,ranges,labels,do_title=True,fsize=FS) :
    nb=128

    sigma_max=0
    i1=find_param(params,name)
    title="$\\sigma($"+params[i1].label+"$)=["
    for i in np.arange(len(fishermat)) :
        covar_full=np.linalg.inv(fishermat[i])
        sigma=np.sqrt(covar_full[i1,i1])
        if sigma>=sigma_max :
            sigma_max=sigma
        x_arr=params[i1].val-4*sigma+8*sigma*np.arange(nb)/(nb-1.)
        p_arr=np.exp(-(x_arr-params[i1].val)**2/(2*sigma**2))
        title+="%.3lf"%sigma
        if i<(len(fishermat)-1) :
            title+=","
#        ax.plot(x_arr,p_arr,color=lc[i],linestyle=ls[i],linewidth=lw[i])
        ax.plot(x_arr,p_arr,color=lc[i],linestyle=ls[i],linewidth=lw[i],label=labels[i])
    title+="]$"
    if do_title :
        ax.set_title(title)
    else :
        print title
    if ranges[name]==None :
        ax.set_xlim([params[i1].val-fact_axis*sigma_max,params[i1].val+fact_axis*sigma_max])
    else :
        ax.set_xlim(ranges[name])
    ax.set_ylim([0,1.02])
    ax.set_xlabel(params[i1].label,fontsize=FS)
    for label in ax.get_yticklabels():
        label.set_fontsize(fsize-2)
    for label in ax.get_xticklabels():
        label.set_fontsize(fsize-2)

def plot_fisher_two(params,name1,name2,fishermat,ax,fc,lw,ls,lc,fact_axis,ranges,fsize=FS) :
    sig0_max=0
    sig1_max=0
    for i in np.arange(len(fishermat)) :
        i1=find_param(params,name1)
        i2=find_param(params,name2)
        covar_full=np.linalg.inv(fishermat[i])
        covar=np.zeros([2,2])
        covar[0,0]=covar_full[i1,i1]
        covar[0,1]=covar_full[i1,i2]
        covar[1,0]=covar_full[i2,i1]
        covar[1,1]=covar_full[i2,i2]
        sig0=np.sqrt(covar[0,0])
        sig1=np.sqrt(covar[1,1])

        if sig0>=sig0_max :
            sig0_max=sig0
        if sig1>=sig1_max :
            sig1_max=sig1

        w,v=np.linalg.eigh(covar)
        angle=180*np.arctan2(v[1,0],v[0,0])/np.pi
        a_1s=np.sqrt(2.3*w[0])
        b_1s=np.sqrt(2.3*w[1])
        a_2s=np.sqrt(6.17*w[0])
        b_2s=np.sqrt(6.17*w[1])

        centre=np.array([params[i1].val,params[i2].val])

        e_1s=Ellipse(xy=centre,width=2*a_1s,height=2*b_1s,angle=angle,
                     facecolor=fc[i],linewidth=lw[i],linestyle=ls[i],edgecolor=lc[i])
#        e_2s=Ellipse(xy=centre,width=2*a_2s,height=2*b_2s,angle=angle,
#                     facecolor=fc[i],linewidth=lw[i],linestyle='dashed',edgecolor=lc[i])

#        ax.add_artist(e_2s)
        ax.add_artist(e_1s)
        if ranges[name1]==None :
            ax.set_xlim([params[i1].val-fact_axis*sig0_max,
                         params[i1].val+fact_axis*sig0_max])
        else :
            ax.set_xlim(ranges[name1])
        if ranges[name2]==None :
            ax.set_ylim([params[i2].val-fact_axis*sig1_max,
                         params[i2].val+fact_axis*sig1_max])
        else :
            ax.set_ylim(ranges[name2])
        ax.set_xlabel(params[i1].label,fontsize=FS)
        ax.set_ylabel(params[i2].label,fontsize=FS)
    for label in ax.get_yticklabels():
        label.set_fontsize(fsize-2)
    for label in ax.get_xticklabels():
        label.set_fontsize(fsize-2)

def plot_fisher_two_standalone(params,name1,name2,fishermat,fc,lw,lc,lab,fact_axis,ranges,fsize=FS) :
    sig0_max=0
    sig1_max=0
    leg_items=[]
    for i in np.arange(len(fishermat)) :
        i1=find_param(params,name1)
        i2=find_param(params,name2)
        covar_full=np.linalg.inv(fishermat[i])
        covar=np.zeros([2,2])
        covar[0,0]=covar_full[i1,i1]
        covar[0,1]=covar_full[i1,i2]
        covar[1,0]=covar_full[i2,i1]
        covar[1,1]=covar_full[i2,i2]
        sig0=np.sqrt(covar[0,0])
        sig1=np.sqrt(covar[1,1])

        if sig0>=sig0_max :
            sig0_max=sig0
        if sig1>=sig1_max :
            sig1_max=sig1

        w,v=np.linalg.eigh(covar)
        angle=180*np.arctan2(v[1,0],v[0,0])/np.pi
        a_1s=np.sqrt(2.3*w[0])
        b_1s=np.sqrt(2.3*w[1])
        a_2s=np.sqrt(6.17*w[0])
        b_2s=np.sqrt(6.17*w[1])

        centre=np.array([params[i1].val,params[i2].val])

        e_1s=Ellipse(xy=centre,width=2*a_1s,height=2*b_1s,angle=angle,
                     facecolor=fc[i],linewidth=lw[i],linestyle="solid",edgecolor=lc[i])
        e_2s=Ellipse(xy=centre,width=2*a_2s,height=2*b_2s,angle=angle,
                     facecolor=fc[i],linewidth=lw[i],linestyle="dashed",edgecolor=lc[i])
        
        ax=plt.gca()
        ax.add_artist(e_2s)
        ax.add_artist(e_1s)
        if ranges[name1]==None :
            ax.set_xlim([params[i1].val-fact_axis*sig0_max,
                         params[i1].val+fact_axis*sig0_max])
        else :
            ax.set_xlim(ranges[name1])
        if ranges[name2]==None :
            ax.set_ylim([params[i2].val-fact_axis*sig1_max,
                         params[i2].val+fact_axis*sig1_max])
        else :
            ax.set_ylim(ranges[name2])
        ax.set_xlabel(params[i1].label,fontsize=FS)
        ax.set_ylabel(params[i2].label,fontsize=FS)
        leg_items.append(plt.Line2D((0,1),(0,0),color=lc[i],linewidth=lw[i]))
    for label in ax.get_yticklabels():
        label.set_fontsize(fsize-2)
    for label in ax.get_xticklabels():
        label.set_fontsize(fsize-2)
    ax.legend(leg_items,lab,loc='upper left',frameon=False)


def plot_fisher_all(params, #Parameters in the FMs
                    fishermat, #FMs to plot
                    fc,lw,ls,lc, #Foreground colors, line widths, line styles and line colours for each FM
                    labels, #Labels for each FM
                    fact_axis, #The x and y axes will be fact_axis x error in each parameter
                    ranges, #Ranges for each parameter
                    fname,
                    do_1D=True,
                    do_titles=True,
                    fsize=FS) : #File to save the plot
    index_plot=np.where(np.array([p.do_plot for p in params]))
    n_params=len(index_plot[0])
    param_plot=params[index_plot]

    fig=plt.figure(figsize=(10,9))
    plt.subplots_adjust(hspace=0,wspace=0)
    for i in np.arange(n_params) : #Plot pdfs and ellipses
        i_col=i
        for j in np.arange(n_params-i)+i :
            i_row=j
            if do_1D :
                iplot=i_col+n_params*i_row+1
            else :
                iplot=i_col+(n_params-1)*(i_row-1)+1

            ax=None
            if do_1D :
                ax=fig.add_subplot(n_params,n_params,iplot)
                if i==j :
                    plot_fisher_single(params,param_plot[i].name,fishermat,
                                       ax,fc,lw,ls,lc,fact_axis,ranges,labels,do_titles,fsize=fsize)
            if i!=j :
                if do_1D :
                    ax=fig.add_subplot(n_params,n_params,iplot)
                else :
                    ax=fig.add_subplot(n_params-1,n_params-1,iplot)
                plot_fisher_two(params,param_plot[i].name,param_plot[j].name,
                                fishermat,ax,fc,lw,ls,lc,fact_axis,ranges,fsize=fsize)

            if ax!=None :
                if i_row!=n_params-1 :
                    ax.get_xaxis().set_visible(False)

                if i_col!=0 :
                    ax.get_yaxis().set_visible(False)

                if i_col==0 and i_row==0 :
                    ax.get_yaxis().set_visible(False)
                
                ax.locator_params(nbins=6)

    if n_params!=2 :
        if n_params>1 : #Add labels in a separate plot
            if do_1D :
                ax=fig.add_subplot(n_params,n_params,2)
            else :
                ax=fig.add_subplot(n_params-1,n_params-1,2)
            ax.set_xlim([-1,1])
            ax.set_ylim([-1,1])
            for i in np.arange(len(labels)) :
                ax.plot([-1,1],[-3,-3],color=lc[i],linestyle=ls[i],
                        linewidth=lw[i],label=labels[i])
            ax.legend(loc='upper left',frameon=False,fontsize=FS)
            ax.axis('off')
        else :
            ax.legend(loc='upper right',frameon=False,fontsize=FS)
#        ax.axis('off')

    if fname!="none" :
        plt.savefig(fname,bbox_inches='tight')

    plt.show()
