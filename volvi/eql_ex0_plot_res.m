function eql_ex0_plot_res(varargin)
    %% *SET-UP*
    inp = varargin{1};
    out = varargin{2};
    set(0,'defaultaxescolororder',[1 0 0;0 0 1]);
    %
    xlm.thv = [0,5];
    [xtk.thv,~] = get_axis_tick(xlm.thv,xlm.thv,inp.vtd/10,1);
    ylm.thv = [-0.020;0.020];
    ytk.thv = -0.020:0.010:0.020;
    %
    xlm.fsv = [1e-1;25];
    xtk.fsv = [1e-1;1e0;1e1;25];
    ylm.fsv = [1e-12;1e0];
    ytk.fsv = [1e-12;1e-9;1e-6;1e0];
    
    [out.vfr,out.fsv(1,:)] = super_fft(mean(diff(out.vtm)),out.thv(1,:),0,[1,2]);
    [out.vfr,out.fsv(2,:)] = super_fft(mean(diff(out.vtm)),out.thv(end,:),0,[1,2]);
    
    %% *PLOT INPUT MOTION*
    data = importdata('thv_volvi_el_strong_bedrock.txt');
    nln.vtm.bdr = data(:,1);
    nln.thv.bdr = data(:,2);
    nln.tha.bdr = [0;diff(nln.thv.bdr)./mean(diff(nln.vtm.bdr))];

    [nln.vfr.bdr,nln.fsv.bdr] = super_fft(mean(diff(nln.vtm.bdr)),nln.thv.bdr,0,[1,2]);
    % time-histories
    fpplot('xpl',{nln.vtm.bdr;out.vtm},'ypl',{nln.thv.bdr;out.thv(end,:)},...
        'vfg','on','pfg',[0,0,14,5.75],'leg',{{'nln','eql'}},'lst',{'-','--'},'lwd',[3,1.5],...
        'xlm',{xlm.thv},'xtk',{xtk.thv},'ylm',{ylm.thv},'ytk',{ytk.thv},...
        'tit',{'G.L. -300m - Strong Motion'},...
        'xlb',{'t [s]'},'ylb',{'v(t) [m/s]'});
    % fourier's spectra
    fpplot('xpl',{nln.vfr.bdr;out.vfr},'ypl',{nln.fsv.bdr;out.fsv(end,:)},...
        'vfg','on','pfg',[0,0,12,12],'leg',{{'nln','eql'}},'lst',{'-','--'},'lwd',[3,1.5],...
        'xlm',{xlm.fsv},'xtk',{xtk.fsv},'ylm',{ylm.fsv},'ytk',{ytk.fsv},...
        'tit',{'G.L. -300m - Strong Motion'},'scl',{'log'},...
        'xlb',{'f [Hz]'},'ylb',{'FSV [m]'});
    
    %% *PLOT OUTPUT MOTION*
    data = importdata('thv_volvi_el_strong_surface.txt');
    nln.vtm.srf = data(:,1);
    nln.thv.srf = data(:,2);
    nln.tha.srf = [0;diff(nln.thv.srf)./mean(diff(nln.vtm.srf))];
    [nln.vfr.srf,nln.fsv.srf] = super_fft(mean(diff(nln.vtm.srf)),nln.thv.srf,0,[1,2]);

    % time-histories
    fpplot('xpl',{nln.vtm.srf;out.vtm},'ypl',{nln.thv.srf;out.thv(1,:)},...
        'vfg','on','pfg',[0,0,14,5.75],'leg',{{'nln','eql'}},'lst',{'-','--'},'lwd',[3,1.5],...
        'xlm',{xlm.thv},'xtk',{xtk.thv},'ylm',{ylm.thv},'ytk',{ytk.thv},...
        'tit',{'G.L. 0m - Strong Motion'},...
        'xlb',{'t [s]'},'ylb',{'v(t) [m/s]'});
    % fourier's spectra
    fpplot('xpl',{nln.vfr.srf;out.vfr},'ypl',{nln.fsv.srf;out.fsv(1,:)},...
        'vfg','on','pfg',[0,0,12,12],'leg',{{'nln','eql'}},'lst',{'-','--'},'lwd',[3,1.5],...
        'xlm',{xlm.fsv},'xtk',{xtk.fsv},'ylm',{ylm.fsv},'ytk',{ytk.fsv},...
        'tit',{'G.L. 0m - Strong Motion'},'scl',{'log'},...
        'xlb',{'f [Hz]'},'ylb',{'FSV [m]'});
    
    %% *PLOT EBHSR MOTION*
    out.bsr = spectral_ratio(out.vtm,out.tha(1,:),out.vtm,out.tha(end,:));
    nln.bsr = spectral_ratio(nln.vtm.srf,nln.tha.srf,nln.vtm.bdr,nln.tha.bdr);
    % fourier's spectra
    fpplot('xpl',{nln.bsr(1,:);out.bsr(1,:)},'ypl',{abs(nln.bsr(2,:));abs(out.bsr(2,:))},...
        'vfg','on','pfg',[0,0,12,12],'leg',{{'nln','eql'}},'lst',{'-','--'},'lwd',[3,1.5],...
        'xlm',{xlm.fsv},'xtk',{xtk.fsv},'ylm',{[1,1e3]},'ytk',{[1,1e1,1e2,1e3]},...
        'tit',{'G.L. 0m/G.L. -300m - Strong Motion'},'scl',{'log'},...
        'xlb',{'f [Hz]'},'ylb',{'BHSR [1]'});
    return
end

