clc;clear;

% makes Phasor diagrams for diff sag types
E=0.4;
V=0.2;

% pre-fault phasors
x1 = [0, 1]*E + [0.5, 0.5];
y1 = [0, 0]*E + [0.5, 0.5];
x2 = [0, -0.5]*E + [0.5, 0.5];
y2 = [0, sqrt(3)/2]*E + [0.5, 0.5];
x3 = [0, -0.5]*E + [0.5, 0.5];
y3 = [0, -sqrt(3)/2]*E + [0.5, 0.5];


va=[E, V, V, E, V, E, V, (2*E/3+V/3)];
vb=[-0.5*E-E*0.5*1j*sqrt(3), -0.5*V-V*0.5*1j*sqrt(3), -0.5*E-E*0.5*1j*sqrt(3), -0.5*E-V*0.5*1j*sqrt(3), -0.5*V-E*0.5*1j*sqrt(3), -0.5*V-V*0.5*1j*sqrt(3), -0.5*V-1j*(V*sqrt(3)/6+E*sqrt(3)/3), -E/3-V/6 - 1j*0.5*sqrt(3)*V];
vc=[-0.5*E+E*0.5*1j*sqrt(3), -0.5*V+V*0.5*1j*sqrt(3), -0.5*E+E*0.5*1j*sqrt(3), -0.5*E+V*0.5*1j*sqrt(3), -0.5*V+E*0.5*1j*sqrt(3), -0.5*V+V*0.5*1j*sqrt(3), -0.5*V+1j*(V*sqrt(3)/6+E*sqrt(3)/3), -E/3-V/6 + 1j*0.5*sqrt(3)*V];
for i=1:8
    fhandle=figure(i); clf; xlim([-1,1]); ylim([-1, 1]); axis off
    set(fhandle, 'Position', [1000, 300, 180, 180])
    annotation('arrow',x1,y1, 'lineStyle', '--')
    annotation('arrow',x2,y2, 'lineStyle', '--')
    annotation('arrow',x3,y3, 'lineStyle', '--')
    
    
    annotation('textarrow',[0, real(va(i))]+ [0.5, 0.5],[0, imag(va(i))]+ [0.5, 0.5], 'color', 'red')
    annotation('textarrow',[0, real(vb(i))]+ [0.5, 0.5],[0, imag(vb(i))]+ [0.5, 0.5], 'color', 'red')
    annotation('textarrow',[0, real(vc(i))]+ [0.5, 0.5],[0, imag(vc(i))]+ [0.5, 0.5], 'color', 'red')
    
    text(real(va(i))*2+0.1, imag(va(i))*2+0.1, 'V_a', 'color', 'red')
    text(real(vb(i))*2-0.5, imag(vb(i))*2, 'V_b', 'color', 'red')
    text(real(vc(i))*2-0.5, imag(vc(i))*2, 'V_c', 'color', 'red')


end
%%
close all