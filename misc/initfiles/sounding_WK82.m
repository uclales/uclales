%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                        %
%  Initial sounding for the squall line described in                     %
%  Weisman and Klemp, 1982, MWR                                          %
%                                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all

g     = 9.81;
cp    = 1005;
th_tr = 343;
z_tr  = 12000;
T_tr  = 213;
th0   = 300;
z_s   = 3000;

z=0:100:18125;
nz=length(z);

% potential temperature th and relative humidity H
for i=1:nz
    if (z(i)<=z_tr)
        th(i)=th0+(th_tr-th0)*(z(i)/z_tr).^(5/4);
        H(i)=1-3/4*(z(i)/z_tr).^(5/4);
    else
        th(i)=th_tr*exp(g/(cp*T_tr)*(z(i)-z_tr));
        H(i)=0.25;
    end
end


% original WK82 wind profile
%u_s = 5.;
%u_s = 15.;
%u_s = 25.;
u_s = 35.;
%u_s = 45.;
u=u_s*tanh(z/z_s);

% set v to 0 m/s
v(1:nz)=0;

fh=figure;
plot(th,z);

fh=figure;
plot(H,z);

fh=figure;
hold on
plot(u,z);
plot(v,z);

z(1)=1000; % 1st line: surface pressure in hPa

% write to file

dlmwrite('bubble_in', [z',th',H',u',v'],'precision','%11.4f','delimiter','','newline', 'unix');


