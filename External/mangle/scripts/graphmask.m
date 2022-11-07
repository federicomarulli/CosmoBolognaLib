function graphmask(infile,outfile,maprange,plottitle,outlines)
% © M E C Swanson 2008
%function for plotting MANGLE polygon files using the Matlab mapping toolbox.  
%arguments: 
%infile=polygon file to be plotted, in 'list' format, i.e.created with poly2poly -ol30 mypolys.pol mypolys.list
%outfile=name of eps file to output, or use 'none' if you want to, e.g., put the resulting figure into a subplot.
%maprange=[lonmin,lonmax,latmin,latmax] are optional latitude and longitude (Dec and RA)
%limits for the mask.  If not provided, default is to plot full sky. 
%title=optional title for the plot
%outlines=whether to draw black outlines around polygons. default is no outlines.  Use outlines='on' to draw outlines.

%check if mapping toolbox is installed:
if ( ~exist('ispolycw') )
     fprintf(2, 'To plot mangle masks in Matlab, you must have the mapping\ntoolbox (version 2.0.3 (R14SP1) or later) installed.\n');
     exitval=1;
     save('matlabexit.temp','exitval');
     exit
end

%process input arguments    
if(nargin<1)
     error('graphmask requires at least one input argument: the name of the input data file to plot.');
end
if(nargin==1)
   outfile='none';
   maprange=[0 360 -90 90];
   plottitle='';
   outlines='';
   lims=0;
end
if(nargin==2)
   maprange=[0 360 -90 90];
   plottitle='';
   outlines='';
   lims=0;
end
if(nargin==3)
   plottitle='';
   outlines='';
   lims=1;
end
if(nargin==4)
   outlines='';
   lims=1;
end
if(nargin>=5)
   lims=1;
end
if(all(maprange==0))
   maprange=[0 360 -90 90];
   lims=0;
end
%check if mapping toolbox is installed:
if (~isnumeric(maprange))
     fprintf(2, 'ERROR: Non-numeric values used in maprange\n');
     exitval=1;
     save('matlabexit.temp','exitval');
     exit
end

lonmin=maprange(1);
lonmax=maprange(2);
latmin=maprange(3);
latmax=maprange(4);

%set spacing for tickmarks
range=max(latmax-latmin, lonmax-lonmin);
avrange=mean([latmax-latmin, lonmax-lonmin]);
if(range<15)
    sp=1;
else if (range<30)
        sp=2;
    else if (range<90)
            sp=5;
        else if (range<150)
                sp=10;
            else if (range<300)
                    sp=20;
                else
                    sp=50;
                end
            end
        end
    end
end
        
%read in files
weightfile=[infile '.weight'];
xymat=load(infile);
wmat=load(weightfile);
fprintf(1,'Done reading files\n');
ra=xymat(1:end,1);
dec=xymat(1:end,2);
weight=wmat(1:end,2);
%set up map axes
axm=axesm('MapProjection', 'hammer','frame','on','FFaceColor','black','origin',180);

if (lims)
    axesm('MapProjection','mercator','frame','on','maplatlimit',[latmin latmax], 'maplonlimit',[lonmin lonmax])
    setm(gca,'ParallelLabel','on','PLabelLocation',sp,'LabelFormat','none','fontsize',10);
    tightmap;
    ax1=gca;
    ax2=axes('Position',get(ax1,'Position'));
    axesm('MapProjection','mercator','frame','on','maplatlimit',[latmin latmax], 'maplonlimit',[lonmin lonmax])
    setm(gca,'FFaceColor','black')
else
   axm=axesm('MapProjection', 'hammer','frame','on','FFaceColor','black','origin',180);
end

%plot polygons in list as patches
if (strcmp(outlines,'on'))
    h=patchesm(dec,ra,'g','Edgecolor','black','Linewidth',0.3);
else
    h=patchesm(dec,ra,'g','Edgecolor','none');
end
%normalize weights
minweight=min(weight);
maxweight=max(weight);
weight=(weight-minweight)./(maxweight-minweight);
%scale weight so that the minimum isn't black unless it's actually zero
if(minweight~=0)
     minscale=.1;
     zeroscale=minweight-minscale./(1-minscale).*(maxweight-minweight);
     if(minscale>0 && zeroscale<=0)
        zeroscale=0;
	minscale=minweight./maxweight;
     end
     weight=(1-minscale).*weight+minscale;
else
     zeroscale=minweight;
end

%set holes (polygons with points wound counterclockwise) to have weight 0
ccw=~ispolycw(dec,ra);
weight(ccw)=0;
%create cell array of colors with each element as grayscale weight
color=[weight weight weight];
cellcolor=num2cell(color,2);
%apply weight colors to patch objects
set(h,{'FaceColor'},cellcolor);

%tweak map display and add labels
if (lims)
    set(gca,'XDir', 'reverse')
    gridm on;
    setm(gca,'MeridianLabel','on','MLabelLocation',sp,'MLineLocation',sp,'PLineLocation',sp,'MLabelParallel','south',...
        'LabelFormat','none', 'MLineLimit',[latmin latmin+.01*avrange],'PLineLimit',[lonmax-.01*avrange, lonmax],...
        'glinestyle','-','gcolor',[.5 .5 .5],'fontsize',10);
    tightmap;
    mlabelzero22pi
    xlabel('\newline Right Ascension')
    ylabel('Declination\newline ')
else
    tightmap;
    xlims=get(gca,'XLim');
    ylims=get(gca,'YLim');
    axis([1.01*xlims 1.01*ylims])
    set(gca,'XDir', 'reverse','XColor',[1 1 1],'YColor',[1 1 1])
    xlabel('-90^{\circ}','Color','k')
    ylabel('360^{\circ}','Rotation',0.0,'Color','k','VerticalAlignment','Middle','HorizontalAlignment','Right')
    ax1=gca;
    ax2=copyobj(ax1,gcf);
    set(ax2,'XAxisLocation','top', 'YAxisLocation','right','Color','none','XColor',[1 1 1],'YColor',[1 1 1]);
    axes(ax2);
    xlabel('+90^{\circ}','Color','k')
    ylabel('0^{\circ}','Rotation',0.0,'Color','k','VerticalAlignment','Middle','HorizontalAlignment','Left')
end
title(plottitle);

%add colorbar
pb=pbaspect;
pos1=get(gca,'Position');
axes('Position',pos1)
axis off
colormap('gray')
caxis([zeroscale maxweight]);
pbaspect(pb);
pos1=get(gca,'Position');
cbar_handle=colorbar('EastOutside');
set(gca,'Position',pos1);

fprintf(1,'Done making mask image\n');
%export image as eps file
if(~strcmp(outfile,'none'))
    set (gcf, 'Color', [1 1 1])
    if(~lims)
        set(gcf,'PaperPosition', [-1.75 .75 14 7])
    end
    print('-depsc','-r600',outfile);
    fprintf(1,'Done writing mask image to %s\n',outfile);
    exit
end 
