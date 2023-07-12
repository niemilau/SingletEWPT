/**** Self energy computations. Will take real parts only ****/
// Some diagrams here are of form d * integral, with d = 4-2eps. So some divergent parts do contribute finite bits.

#include "renormalization.h"

double Renormalization::CalcSelfEnergyH1(double k) {

	// Part proportional to epsilon
	Integral epsPart = (
		-1.0*(ct*ct*g0sq*A0(MW)) - (ct*ct*g0sq*(MZ*MZ)*A0(MZ))/(2.*(MW*MW)) - 2*(ct*ct)*g0sq*(MW*MW)*B0(k,MW,MW) - 
		(ct*ct*g0sq*(MZ*MZ*MZ*MZ)*B0(k,MZ,MZ))/(MW*MW) + (9*(ct*ct)*sqrt(g0sq)*A0(MW)*sqrt(g0sq))/4. + 
		(3*a2*(ct*ct)*(MW*MW)*A0(MW)*sqrt(g0sq))/(sqrt(g0sq)*(Mh1*Mh1)) + (sqrt(g0sq)*(st*st)*A0(MW)*sqrt(g0sq))/2. + 
		(sqrt(g0sq)*(Mh1*Mh1)*(st*st)*A0(MW)*sqrt(g0sq))/(Mh2*Mh2) - 
		(2*a2*(MW*MW)*(st*st)*A0(MW)*sqrt(g0sq))/(sqrt(g0sq)*(Mh2*Mh2)) - 
		(4*b3*ct*MW*(st*st*st)*A0(MW)*sqrt(g0sq))/(Mh1*Mh1) + 
		(3*a2*(ct*ct)*(MZ*MZ)*A0(MZ)*sqrt(g0sq))/(2.*sqrt(g0sq)*(Mh1*Mh1)) + 
		(9*(ct*ct)*sqrt(g0sq)*(MZ*MZ)*A0(MZ)*sqrt(g0sq))/(8.*(MW*MW)) - 
		(a2*(MZ*MZ)*(st*st)*A0(MZ)*sqrt(g0sq))/(sqrt(g0sq)*(Mh2*Mh2)) + 
		(sqrt(g0sq)*(MZ*MZ)*(st*st)*A0(MZ)*sqrt(g0sq))/(4.*(MW*MW)) + 
		(sqrt(g0sq)*(Mh1*Mh1)*(MZ*MZ)*(st*st)*A0(MZ)*sqrt(g0sq))/(2.*(Mh2*Mh2)*(MW*MW)) - 
		(2*b3*ct*(MZ*MZ)*(st*st*st)*A0(MZ)*sqrt(g0sq))/(Mh1*Mh1*MW) + 
		(sqrt(g0sq)*(st*st)*A0(MW)*cos(2*theta)*sqrt(g0sq))/2. + 
		(sqrt(g0sq)*(Mh1*Mh1)*(st*st)*A0(MW)*cos(2*theta)*sqrt(g0sq))/(Mh2*Mh2) - 
		(6*a2*(MW*MW)*(st*st)*A0(MW)*cos(2*theta)*sqrt(g0sq))/(sqrt(g0sq)*(Mh2*Mh2)) - 
		(3*a2*(MZ*MZ)*(st*st)*A0(MZ)*cos(2*theta)*sqrt(g0sq))/(sqrt(g0sq)*(Mh2*Mh2)) + 
		(sqrt(g0sq)*(MZ*MZ)*(st*st)*A0(MZ)*cos(2*theta)*sqrt(g0sq))/(4.*(MW*MW)) + 
		(sqrt(g0sq)*(Mh1*Mh1)*(MZ*MZ)*(st*st)*A0(MZ)*cos(2*theta)*sqrt(g0sq))/(2.*(Mh2*Mh2)*(MW*MW)) + 
		(3*ct*sqrt(g0sq)*A0(MW)*cos(3*theta)*sqrt(g0sq))/4. - 
		(3*a2*ct*(MW*MW)*A0(MW)*cos(3*theta)*sqrt(g0sq))/(sqrt(g0sq)*(Mh1*Mh1)) - 
		(3*a2*ct*(MZ*MZ)*A0(MZ)*cos(3*theta)*sqrt(g0sq))/(2.*sqrt(g0sq)*(Mh1*Mh1)) + 
		(3*ct*sqrt(g0sq)*(MZ*MZ)*A0(MZ)*cos(3*theta)*sqrt(g0sq))/(8.*(MW*MW)) + 
		(2*b3*MW*(st*st)*A0(MW)*sin(2*theta)*sqrt(g0sq))/(Mh2*Mh2) + 
		(b3*(MZ*MZ)*(st*st)*A0(MZ)*sin(2*theta)*sqrt(g0sq))/(Mh2*Mh2*MW)
	);

	Integral others = (
        (3*(ct*ct*ct*ct)*g0sq*(Mh1*Mh1)*A0(Mh1))/(16.*(MW*MW)) + (3*(ct*ct*ct*ct)*g0sq*(Mh2*Mh2)*A0(Mh1))/(16.*(MW*MW)) + 
        3*b4*(st*st*st*st)*A0(Mh1) + (a2*A0(Mh2))/8. + (3*b4*A0(Mh2))/8. + (3*g0sq*(Mh1*Mh1)*A0(Mh2))/(128.*(MW*MW)) + 
        (3*g0sq*(Mh2*Mh2)*A0(Mh2))/(128.*(MW*MW)) + 2*(ct*ct)*g0sq*A0(MW) + (ct*ct*g0sq*(Mh1*Mh1)*A0(MW))/(8.*(MW*MW)) + 
        (ct*ct*g0sq*(Mh2*Mh2)*A0(MW))/(8.*(MW*MW)) + a2*(st*st)*A0(MW) + (ct*ct*g0sq*(Mh1*Mh1)*A0(MZ))/(16.*(MW*MW)) + 
        (ct*ct*g0sq*(Mh2*Mh2)*A0(MZ))/(16.*(MW*MW)) + (ct*ct*g0sq*(MZ*MZ)*A0(MZ))/(MW*MW) + (a2*(st*st)*A0(MZ))/2. + 
        (27*a2*(ct*ct)*(Mh1*Mh1)*B0(k,Mh1,Mh1))/16. + (81*(ct*ct)*g0sq*(Mh1*Mh1*Mh1*Mh1)*B0(k,Mh1,Mh1))/(128.*(MW*MW)) + 
        (9*(a2*a2)*(ct*ct)*(MW*MW)*B0(k,Mh1,Mh1))/(8.*g0sq) - 
        (9*b3*ct*sqrt(g0sq)*(Mh1*Mh1)*(st*st*st)*B0(k,Mh1,Mh1))/(4.*MW) - 
        (3*a2*b3*ct*MW*(st*st*st)*B0(k,Mh1,Mh1))/sqrt(g0sq) + 2*(b3*b3)*(st*st*st*st*st*st)*B0(k,Mh1,Mh1) - 
        a2*(Mh1*Mh1)*(st*st)*B0(k,Mh1,Mh2) - (a2*(Mh2*Mh2)*(st*st)*B0(k,Mh1,Mh2))/2. + 
        (g0sq*(Mh1*Mh1*Mh1*Mh1)*(st*st)*B0(k,Mh1,Mh2))/(4.*(MW*MW)) + 
        (g0sq*(Mh1*Mh1)*(Mh2*Mh2)*(st*st)*B0(k,Mh1,Mh2))/(4.*(MW*MW)) + 
        (g0sq*(Mh2*Mh2*Mh2*Mh2)*(st*st)*B0(k,Mh1,Mh2))/(16.*(MW*MW)) + (a2*a2*(MW*MW)*(st*st)*B0(k,Mh1,Mh2))/g0sq - 
        (a2*(ct*ct)*(Mh1*Mh1)*B0(k,Mh2,Mh2))/4. - (a2*(ct*ct)*(Mh2*Mh2)*B0(k,Mh2,Mh2))/2. + 
        (ct*ct*g0sq*(Mh1*Mh1*Mh1*Mh1)*B0(k,Mh2,Mh2))/(32.*(MW*MW)) + 
        (ct*ct*g0sq*(Mh1*Mh1)*(Mh2*Mh2)*B0(k,Mh2,Mh2))/(8.*(MW*MW)) + 
        (ct*ct*g0sq*(Mh2*Mh2*Mh2*Mh2)*B0(k,Mh2,Mh2))/(8.*(MW*MW)) + (a2*a2*(ct*ct)*(MW*MW)*B0(k,Mh2,Mh2))/(2.*g0sq) + 
        (ct*ct*g0sq*(Mh1*Mh1*Mh1*Mh1)*B0(k,MW,MW))/(4.*(MW*MW)) + (7*(ct*ct)*g0sq*(MW*MW)*B0(k,MW,MW))/2. + 
        (ct*ct*g0sq*(Mh1*Mh1*Mh1*Mh1)*B0(k,MZ,MZ))/(8.*(MW*MW)) + (7*(ct*ct)*g0sq*(MZ*MZ*MZ*MZ)*B0(k,MZ,MZ))/(4.*(MW*MW)) + 
        (3*(ct*ct*ct*ct)*g0sq*(Mh1*Mh1)*A0(Mh1)*cos(2*theta))/(16.*(MW*MW)) - 
        (3*(ct*ct*ct*ct)*g0sq*(Mh2*Mh2)*A0(Mh1)*cos(2*theta))/(16.*(MW*MW)) + 
        (3*g0sq*(Mh1*Mh1)*A0(Mh2)*cos(2*theta))/(256.*(MW*MW)) - (3*g0sq*(Mh2*Mh2)*A0(Mh2)*cos(2*theta))/(256.*(MW*MW)) + 
        (ct*ct*g0sq*(Mh1*Mh1)*A0(MW)*cos(2*theta))/(8.*(MW*MW)) - (ct*ct*g0sq*(Mh2*Mh2)*A0(MW)*cos(2*theta))/(8.*(MW*MW)) + 
        (ct*ct*g0sq*(Mh1*Mh1)*A0(MZ)*cos(2*theta))/(16.*(MW*MW)) - (ct*ct*g0sq*(Mh2*Mh2)*A0(MZ)*cos(2*theta))/(16.*(MW*MW)) - 
        4*a2*(Mh1*Mh1)*(st*st)*B0(k,Mh1,Mh2)*cos(2*theta) - 2*a2*(Mh2*Mh2)*(st*st)*B0(k,Mh1,Mh2)*cos(2*theta) + 
        (g0sq*(Mh1*Mh1*Mh1*Mh1)*(st*st)*B0(k,Mh1,Mh2)*cos(2*theta))/(2.*(MW*MW)) + 
        (g0sq*(Mh1*Mh1)*(Mh2*Mh2)*(st*st)*B0(k,Mh1,Mh2)*cos(2*theta))/(2.*(MW*MW)) + 
        (g0sq*(Mh2*Mh2*Mh2*Mh2)*(st*st)*B0(k,Mh1,Mh2)*cos(2*theta))/(8.*(MW*MW)) + 
        (6*(a2*a2)*(MW*MW)*(st*st)*B0(k,Mh1,Mh2)*cos(2*theta))/g0sq + a2*(ct*ct)*(Mh1*Mh1)*B0(k,Mh2,Mh2)*cos(2*theta) + 
        2*a2*(ct*ct)*(Mh2*Mh2)*B0(k,Mh2,Mh2)*cos(2*theta) - 
        (ct*ct*g0sq*(Mh1*Mh1*Mh1*Mh1)*B0(k,Mh2,Mh2)*cos(2*theta))/(16.*(MW*MW)) - 
        (ct*ct*g0sq*(Mh1*Mh1)*(Mh2*Mh2)*B0(k,Mh2,Mh2)*cos(2*theta))/(4.*(MW*MW)) - 
        (ct*ct*g0sq*(Mh2*Mh2*Mh2*Mh2)*B0(k,Mh2,Mh2)*cos(2*theta))/(4.*(MW*MW)) - 
        (3*(a2*a2)*(ct*ct)*(MW*MW)*B0(k,Mh2,Mh2)*cos(2*theta))/g0sq - 
        3*a2*(Mh1*Mh1)*(st*st)*B0(k,Mh1,Mh2)*(cos(2*theta)*cos(2*theta)) - 
        (3*a2*(Mh2*Mh2)*(st*st)*B0(k,Mh1,Mh2)*(cos(2*theta)*cos(2*theta)))/2. + 
        (g0sq*(Mh1*Mh1*Mh1*Mh1)*(st*st)*B0(k,Mh1,Mh2)*(cos(2*theta)*cos(2*theta)))/(4.*(MW*MW)) + 
        (g0sq*(Mh1*Mh1)*(Mh2*Mh2)*(st*st)*B0(k,Mh1,Mh2)*(cos(2*theta)*cos(2*theta)))/(4.*(MW*MW)) + 
        (g0sq*(Mh2*Mh2*Mh2*Mh2)*(st*st)*B0(k,Mh1,Mh2)*(cos(2*theta)*cos(2*theta)))/(16.*(MW*MW)) + 
        (9*(a2*a2)*(MW*MW)*(st*st)*B0(k,Mh1,Mh2)*(cos(2*theta)*cos(2*theta)))/g0sq - 
        (3*a2*(ct*ct)*(Mh1*Mh1)*B0(k,Mh2,Mh2)*(cos(2*theta)*cos(2*theta)))/4. - 
        (3*a2*(ct*ct)*(Mh2*Mh2)*B0(k,Mh2,Mh2)*(cos(2*theta)*cos(2*theta)))/2. + 
        (ct*ct*g0sq*(Mh1*Mh1*Mh1*Mh1)*B0(k,Mh2,Mh2)*(cos(2*theta)*cos(2*theta)))/(32.*(MW*MW)) + 
        (ct*ct*g0sq*(Mh1*Mh1)*(Mh2*Mh2)*B0(k,Mh2,Mh2)*(cos(2*theta)*cos(2*theta)))/(8.*(MW*MW)) + 
        (ct*ct*g0sq*(Mh2*Mh2*Mh2*Mh2)*B0(k,Mh2,Mh2)*(cos(2*theta)*cos(2*theta)))/(8.*(MW*MW)) + 
        (9*(a2*a2)*(ct*ct)*(MW*MW)*B0(k,Mh2,Mh2)*(cos(2*theta)*cos(2*theta)))/(2.*g0sq) - 
        (9*a2*ct*(Mh1*Mh1)*B0(k,Mh1,Mh1)*cos(3*theta))/8. + 
        (27*ct*g0sq*(Mh1*Mh1*Mh1*Mh1)*B0(k,Mh1,Mh1)*cos(3*theta))/(64.*(MW*MW)) - 
        (9*(a2*a2)*ct*(MW*MW)*B0(k,Mh1,Mh1)*cos(3*theta))/(4.*g0sq) - 
        (3*b3*sqrt(g0sq)*(Mh1*Mh1)*(st*st*st)*B0(k,Mh1,Mh1)*cos(3*theta))/(4.*MW) + 
        (3*a2*b3*MW*(st*st*st)*B0(k,Mh1,Mh1)*cos(3*theta))/sqrt(g0sq) - 
        (9*a2*(Mh1*Mh1)*B0(k,Mh1,Mh1)*(cos(3*theta)*cos(3*theta)))/16. + 
        (9*g0sq*(Mh1*Mh1*Mh1*Mh1)*B0(k,Mh1,Mh1)*(cos(3*theta)*cos(3*theta)))/(128.*(MW*MW)) + 
        (9*(a2*a2)*(MW*MW)*B0(k,Mh1,Mh1)*(cos(3*theta)*cos(3*theta)))/(8.*g0sq) + (3*a2*A0(Mh2)*cos(4*theta))/8. - 
        (3*b4*A0(Mh2)*cos(4*theta))/8. - (3*g0sq*(Mh1*Mh1)*A0(Mh2)*cos(4*theta))/(128.*(MW*MW)) - 
        (3*g0sq*(Mh2*Mh2)*A0(Mh2)*cos(4*theta))/(128.*(MW*MW)) - (3*g0sq*(Mh1*Mh1)*A0(Mh2)*cos(6*theta))/(256.*(MW*MW)) + 
        (3*g0sq*(Mh2*Mh2)*A0(Mh2)*cos(6*theta))/(256.*(MW*MW)) - (3*(ct*ct)*g0sq*(Mt*Mt)*K(k,Mt,Mt,1,1,Mt*Mt))/(MW*MW) - 
        (ct*ct*g0sq*K(k,MW,MW,1,4,4*(k*k)))/2. - (ct*ct*g0sq*(MZ*MZ)*K(k,MZ,MZ,1,4,4*(k*k)))/(4.*(MW*MW)) + 
        (b3*sqrt(g0sq)*(Mh1*Mh1)*(st*st)*B0(k,Mh1,Mh2)*sin(2*theta))/MW + 
        (b3*sqrt(g0sq)*(Mh2*Mh2)*(st*st)*B0(k,Mh1,Mh2)*sin(2*theta))/(2.*MW) - 
        (2*a2*b3*MW*(st*st)*B0(k,Mh1,Mh2)*sin(2*theta))/sqrt(g0sq) - 
        (b3*(ct*ct)*sqrt(g0sq)*(Mh1*Mh1)*B0(k,Mh2,Mh2)*sin(2*theta))/(4.*MW) - 
        (b3*(ct*ct)*sqrt(g0sq)*(Mh2*Mh2)*B0(k,Mh2,Mh2)*sin(2*theta))/(2.*MW) + 
        (a2*b3*(ct*ct)*MW*B0(k,Mh2,Mh2)*sin(2*theta))/sqrt(g0sq) + 
        (b3*sqrt(g0sq)*(Mh1*Mh1)*(st*st)*B0(k,Mh1,Mh2)*cos(2*theta)*sin(2*theta))/MW + 
        (b3*sqrt(g0sq)*(Mh2*Mh2)*(st*st)*B0(k,Mh1,Mh2)*cos(2*theta)*sin(2*theta))/(2.*MW) - 
        (6*a2*b3*MW*(st*st)*B0(k,Mh1,Mh2)*cos(2*theta)*sin(2*theta))/sqrt(g0sq) + 
        (b3*(ct*ct)*sqrt(g0sq)*(Mh1*Mh1)*B0(k,Mh2,Mh2)*cos(2*theta)*sin(2*theta))/(4.*MW) + 
        (b3*(ct*ct)*sqrt(g0sq)*(Mh2*Mh2)*B0(k,Mh2,Mh2)*cos(2*theta)*sin(2*theta))/(2.*MW) - 
        (3*a2*b3*(ct*ct)*MW*B0(k,Mh2,Mh2)*cos(2*theta)*sin(2*theta))/sqrt(g0sq) + 
        (3*a2*A0(Mh1)*(sin(2*theta)*sin(2*theta)))/4. + b3*b3*(st*st)*B0(k,Mh1,Mh2)*(sin(2*theta)*sin(2*theta)) + 
        (b3*b3*(ct*ct)*B0(k,Mh2,Mh2)*(sin(2*theta)*sin(2*theta)))/2. - (27*(ct*ct)*sqrt(g0sq)*A0(MW)*sqrt(g0sq))/8. - 
        (9*a2*(ct*ct)*(MW*MW)*A0(MW)*sqrt(g0sq))/(2.*sqrt(g0sq)*(Mh1*Mh1)) - (3*sqrt(g0sq)*(st*st)*A0(MW)*sqrt(g0sq))/4. - 
        (3*sqrt(g0sq)*(Mh1*Mh1)*(st*st)*A0(MW)*sqrt(g0sq))/(2.*(Mh2*Mh2)) + 
        (3*a2*(MW*MW)*(st*st)*A0(MW)*sqrt(g0sq))/(sqrt(g0sq)*(Mh2*Mh2)) + 
        (6*b3*ct*MW*(st*st*st)*A0(MW)*sqrt(g0sq))/(Mh1*Mh1) - 
        (9*a2*(ct*ct)*(MZ*MZ)*A0(MZ)*sqrt(g0sq))/(4.*sqrt(g0sq)*(Mh1*Mh1)) - 
        (27*(ct*ct)*sqrt(g0sq)*(MZ*MZ)*A0(MZ)*sqrt(g0sq))/(16.*(MW*MW)) + 
        (3*a2*(MZ*MZ)*(st*st)*A0(MZ)*sqrt(g0sq))/(2.*sqrt(g0sq)*(Mh2*Mh2)) - 
        (3*sqrt(g0sq)*(MZ*MZ)*(st*st)*A0(MZ)*sqrt(g0sq))/(8.*(MW*MW)) - 
        (3*sqrt(g0sq)*(Mh1*Mh1)*(MZ*MZ)*(st*st)*A0(MZ)*sqrt(g0sq))/(4.*(Mh2*Mh2)*(MW*MW)) + 
        (3*b3*ct*(MZ*MZ)*(st*st*st)*A0(MZ)*sqrt(g0sq))/(Mh1*Mh1*MW) - 
        (3*sqrt(g0sq)*(st*st)*A0(MW)*cos(2*theta)*sqrt(g0sq))/4. - 
        (3*sqrt(g0sq)*(Mh1*Mh1)*(st*st)*A0(MW)*cos(2*theta)*sqrt(g0sq))/(2.*(Mh2*Mh2)) + 
        (9*a2*(MW*MW)*(st*st)*A0(MW)*cos(2*theta)*sqrt(g0sq))/(sqrt(g0sq)*(Mh2*Mh2)) + 
        (9*a2*(MZ*MZ)*(st*st)*A0(MZ)*cos(2*theta)*sqrt(g0sq))/(2.*sqrt(g0sq)*(Mh2*Mh2)) - 
        (3*sqrt(g0sq)*(MZ*MZ)*(st*st)*A0(MZ)*cos(2*theta)*sqrt(g0sq))/(8.*(MW*MW)) - 
        (3*sqrt(g0sq)*(Mh1*Mh1)*(MZ*MZ)*(st*st)*A0(MZ)*cos(2*theta)*sqrt(g0sq))/(4.*(Mh2*Mh2)*(MW*MW)) - 
        (9*ct*sqrt(g0sq)*A0(MW)*cos(3*theta)*sqrt(g0sq))/8. + 
        (9*a2*ct*(MW*MW)*A0(MW)*cos(3*theta)*sqrt(g0sq))/(2.*sqrt(g0sq)*(Mh1*Mh1)) + 
        (9*a2*ct*(MZ*MZ)*A0(MZ)*cos(3*theta)*sqrt(g0sq))/(4.*sqrt(g0sq)*(Mh1*Mh1)) - 
        (9*ct*sqrt(g0sq)*(MZ*MZ)*A0(MZ)*cos(3*theta)*sqrt(g0sq))/(16.*(MW*MW)) - 
        (3*b3*MW*(st*st)*A0(MW)*sin(2*theta)*sqrt(g0sq))/(Mh2*Mh2) - 
        (3*b3*(MZ*MZ)*(st*st)*A0(MZ)*sin(2*theta)*sqrt(g0sq))/(2.*(Mh2*Mh2)*MW) - 
        (27*a2*(ct*ct)*A0(Mh1)*(sqrt(g0sq)*sqrt(g0sq)))/(16.*g0sq) - 
        (81*(ct*ct)*(Mh1*Mh1)*A0(Mh1)*(sqrt(g0sq)*sqrt(g0sq)))/(128.*(MW*MW)) - 
        (9*(a2*a2)*(ct*ct)*(MW*MW)*A0(Mh1)*(sqrt(g0sq)*sqrt(g0sq)))/(8.*(g0sq*g0sq)*(Mh1*Mh1)) + 
        (a2*(st*st)*A0(Mh1)*(sqrt(g0sq)*sqrt(g0sq)))/(4.*g0sq) + 
        (a2*(Mh1*Mh1)*(st*st)*A0(Mh1)*(sqrt(g0sq)*sqrt(g0sq)))/(2.*g0sq*(Mh2*Mh2)) - 
        (Mh1*Mh1*(st*st)*A0(Mh1)*(sqrt(g0sq)*sqrt(g0sq)))/(8.*(MW*MW)) - 
        (Mh1*Mh1*Mh1*Mh1*(st*st)*A0(Mh1)*(sqrt(g0sq)*sqrt(g0sq)))/(8.*(Mh2*Mh2)*(MW*MW)) - 
        (Mh2*Mh2*(st*st)*A0(Mh1)*(sqrt(g0sq)*sqrt(g0sq)))/(32.*(MW*MW)) - 
        (a2*a2*(MW*MW)*(st*st)*A0(Mh1)*(sqrt(g0sq)*sqrt(g0sq)))/(2.*(g0sq*g0sq)*(Mh2*Mh2)) + 
        (9*b3*ct*(st*st*st)*A0(Mh1)*(sqrt(g0sq)*sqrt(g0sq)))/(4.*sqrt(g0sq)*MW) + 
        (3*a2*b3*ct*MW*(st*st*st)*A0(Mh1)*(sqrt(g0sq)*sqrt(g0sq)))/(g0sq*sqrt(g0sq)*(Mh1*Mh1)) - 
        (2*(b3*b3)*(st*st*st*st*st*st)*A0(Mh1)*(sqrt(g0sq)*sqrt(g0sq)))/(g0sq*(Mh1*Mh1)) + 
        (3*a2*(ct*ct)*A0(Mh2)*(sqrt(g0sq)*sqrt(g0sq)))/(8.*g0sq) - 
        (3*a2*(ct*ct)*(Mh2*Mh2)*A0(Mh2)*(sqrt(g0sq)*sqrt(g0sq)))/(8.*g0sq*(Mh1*Mh1)) - 
        (9*(ct*ct)*(Mh1*Mh1)*A0(Mh2)*(sqrt(g0sq)*sqrt(g0sq)))/(64.*(MW*MW)) - 
        (9*(ct*ct)*(Mh2*Mh2)*A0(Mh2)*(sqrt(g0sq)*sqrt(g0sq)))/(32.*(MW*MW)) + 
        (3*(a2*a2)*(ct*ct)*(MW*MW)*A0(Mh2)*(sqrt(g0sq)*sqrt(g0sq)))/(4.*(g0sq*g0sq)*(Mh1*Mh1)) - 
        (3*b3*ct*st*A0(Mh2)*(sqrt(g0sq)*sqrt(g0sq)))/(16.*sqrt(g0sq)*MW) - 
        (3*b3*ct*(Mh1*Mh1)*st*A0(Mh2)*(sqrt(g0sq)*sqrt(g0sq)))/(8.*sqrt(g0sq)*(Mh2*Mh2)*MW) + 
        (3*a2*b3*ct*MW*st*A0(Mh2)*(sqrt(g0sq)*sqrt(g0sq)))/(4.*g0sq*sqrt(g0sq)*(Mh2*Mh2)) - 
        (3*a2*(Mh1*Mh1)*(st*st)*A0(Mh2)*(sqrt(g0sq)*sqrt(g0sq)))/(4.*g0sq*(Mh2*Mh2)) - 
        (3*(Mh1*Mh1)*(st*st)*A0(Mh2)*(sqrt(g0sq)*sqrt(g0sq)))/(16.*(MW*MW)) - 
        (3*(Mh2*Mh2)*(st*st)*A0(Mh2)*(sqrt(g0sq)*sqrt(g0sq)))/(32.*(MW*MW)) + 
        (3*(a2*a2)*(MW*MW)*(st*st)*A0(Mh2)*(sqrt(g0sq)*sqrt(g0sq)))/(2.*(g0sq*g0sq)*(Mh2*Mh2)) + 
        (b3*ct*(st*st*st)*A0(Mh2)*(sqrt(g0sq)*sqrt(g0sq)))/(4.*sqrt(g0sq)*MW) + 
        (b3*ct*(Mh2*Mh2)*(st*st*st)*A0(Mh2)*(sqrt(g0sq)*sqrt(g0sq)))/(2.*sqrt(g0sq)*(Mh1*Mh1)*MW) - 
        (a2*b3*ct*MW*(st*st*st)*A0(Mh2)*(sqrt(g0sq)*sqrt(g0sq)))/(g0sq*sqrt(g0sq)*(Mh1*Mh1)) + 
        (9*a2*(ct*ct)*(Mt*Mt)*A0(Mt)*(sqrt(g0sq)*sqrt(g0sq)))/(g0sq*(Mh1*Mh1)) + 
        (27*(ct*ct)*(Mt*Mt)*A0(Mt)*(sqrt(g0sq)*sqrt(g0sq)))/(4.*(MW*MW)) - 
        (6*a2*(Mt*Mt)*(st*st)*A0(Mt)*(sqrt(g0sq)*sqrt(g0sq)))/(g0sq*(Mh2*Mh2)) + 
        (3*(Mt*Mt)*(st*st)*A0(Mt)*(sqrt(g0sq)*sqrt(g0sq)))/(2.*(MW*MW)) + 
        (3*(Mh1*Mh1)*(Mt*Mt)*(st*st)*A0(Mt)*(sqrt(g0sq)*sqrt(g0sq)))/(Mh2*Mh2*(MW*MW)) - 
        (12*b3*ct*(Mt*Mt)*(st*st*st)*A0(Mt)*(sqrt(g0sq)*sqrt(g0sq)))/(sqrt(g0sq)*(Mh1*Mh1)*MW) - 
        (3*a2*(ct*ct)*A0(MW)*(sqrt(g0sq)*sqrt(g0sq)))/(4.*g0sq) - 
        (9*(ct*ct)*(Mh1*Mh1)*A0(MW)*(sqrt(g0sq)*sqrt(g0sq)))/(16.*(MW*MW)) + 
        (a2*(st*st)*A0(MW)*(sqrt(g0sq)*sqrt(g0sq)))/(2.*g0sq) - 
        (Mh1*Mh1*(st*st)*A0(MW)*(sqrt(g0sq)*sqrt(g0sq)))/(4.*(MW*MW)) - 
        (Mh2*Mh2*(st*st)*A0(MW)*(sqrt(g0sq)*sqrt(g0sq)))/(8.*(MW*MW)) + 
        (b3*ct*(st*st*st)*A0(MW)*(sqrt(g0sq)*sqrt(g0sq)))/(sqrt(g0sq)*MW) - 
        (3*a2*(ct*ct)*A0(MZ)*(sqrt(g0sq)*sqrt(g0sq)))/(8.*g0sq) - 
        (9*(ct*ct)*(Mh1*Mh1)*A0(MZ)*(sqrt(g0sq)*sqrt(g0sq)))/(32.*(MW*MW)) + 
        (a2*(st*st)*A0(MZ)*(sqrt(g0sq)*sqrt(g0sq)))/(4.*g0sq) - 
        (Mh1*Mh1*(st*st)*A0(MZ)*(sqrt(g0sq)*sqrt(g0sq)))/(8.*(MW*MW)) - 
        (Mh2*Mh2*(st*st)*A0(MZ)*(sqrt(g0sq)*sqrt(g0sq)))/(16.*(MW*MW)) + 
        (b3*ct*(st*st*st)*A0(MZ)*(sqrt(g0sq)*sqrt(g0sq)))/(2.*sqrt(g0sq)*MW) + 
        (a2*(st*st)*A0(Mh1)*cos(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/g0sq + 
        (2*a2*(Mh1*Mh1)*(st*st)*A0(Mh1)*cos(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(g0sq*(Mh2*Mh2)) - 
        (Mh1*Mh1*(st*st)*A0(Mh1)*cos(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(4.*(MW*MW)) - 
        (Mh1*Mh1*Mh1*Mh1*(st*st)*A0(Mh1)*cos(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(4.*(Mh2*Mh2)*(MW*MW)) - 
        (Mh2*Mh2*(st*st)*A0(Mh1)*cos(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(16.*(MW*MW)) - 
        (3*(a2*a2)*(MW*MW)*(st*st)*A0(Mh1)*cos(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(g0sq*g0sq*(Mh2*Mh2)) - 
        (3*a2*(ct*ct)*A0(Mh2)*cos(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(2.*g0sq) + 
        (3*a2*(ct*ct)*(Mh2*Mh2)*A0(Mh2)*cos(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(8.*g0sq*(Mh1*Mh1)) + 
        (9*(ct*ct)*(Mh1*Mh1)*A0(Mh2)*cos(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(64.*(MW*MW)) + 
        (9*(ct*ct)*(Mh2*Mh2)*A0(Mh2)*cos(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(32.*(MW*MW)) - 
        (9*(a2*a2)*(ct*ct)*(MW*MW)*A0(Mh2)*cos(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(4.*(g0sq*g0sq)*(Mh1*Mh1)) - 
        (3*b3*ct*st*A0(Mh2)*cos(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(16.*sqrt(g0sq)*MW) - 
        (3*b3*ct*(Mh1*Mh1)*st*A0(Mh2)*cos(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(8.*sqrt(g0sq)*(Mh2*Mh2)*MW) + 
        (9*a2*b3*ct*MW*st*A0(Mh2)*cos(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(4.*g0sq*sqrt(g0sq)*(Mh2*Mh2)) - 
        (3*a2*(Mh1*Mh1)*(st*st)*A0(Mh2)*cos(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(2.*g0sq*(Mh2*Mh2)) + 
        (6*(a2*a2)*(MW*MW)*(st*st)*A0(Mh2)*cos(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(g0sq*g0sq*(Mh2*Mh2)) - 
        (b3*ct*(st*st*st)*A0(Mh2)*cos(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(4.*sqrt(g0sq)*MW) - 
        (b3*ct*(Mh2*Mh2)*(st*st*st)*A0(Mh2)*cos(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(2.*sqrt(g0sq)*(Mh1*Mh1)*MW) + 
        (3*a2*b3*ct*MW*(st*st*st)*A0(Mh2)*cos(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(g0sq*sqrt(g0sq)*(Mh1*Mh1)) - 
        (18*a2*(Mt*Mt)*(st*st)*A0(Mt)*cos(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(g0sq*(Mh2*Mh2)) + 
        (3*(Mt*Mt)*(st*st)*A0(Mt)*cos(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(2.*(MW*MW)) + 
        (3*(Mh1*Mh1)*(Mt*Mt)*(st*st)*A0(Mt)*cos(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(Mh2*Mh2*(MW*MW)) + 
        (3*a2*(st*st)*A0(MW)*cos(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(2.*g0sq) - 
        (Mh1*Mh1*(st*st)*A0(MW)*cos(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(4.*(MW*MW)) - 
        (Mh2*Mh2*(st*st)*A0(MW)*cos(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(8.*(MW*MW)) + 
        (3*a2*(st*st)*A0(MZ)*cos(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(4.*g0sq) - 
        (Mh1*Mh1*(st*st)*A0(MZ)*cos(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(8.*(MW*MW)) - 
        (Mh2*Mh2*(st*st)*A0(MZ)*cos(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(16.*(MW*MW)) + 
        (3*a2*(st*st)*A0(Mh1)*(cos(2*theta)*cos(2*theta))*(sqrt(g0sq)*sqrt(g0sq)))/(4.*g0sq) + 
        (3*a2*(Mh1*Mh1)*(st*st)*A0(Mh1)*(cos(2*theta)*cos(2*theta))*(sqrt(g0sq)*sqrt(g0sq)))/(2.*g0sq*(Mh2*Mh2)) - 
        (Mh1*Mh1*(st*st)*A0(Mh1)*(cos(2*theta)*cos(2*theta))*(sqrt(g0sq)*sqrt(g0sq)))/(8.*(MW*MW)) - 
        (Mh1*Mh1*Mh1*Mh1*(st*st)*A0(Mh1)*(cos(2*theta)*cos(2*theta))*(sqrt(g0sq)*sqrt(g0sq)))/(8.*(Mh2*Mh2)*(MW*MW)) - 
        (Mh2*Mh2*(st*st)*A0(Mh1)*(cos(2*theta)*cos(2*theta))*(sqrt(g0sq)*sqrt(g0sq)))/(32.*(MW*MW)) - 
        (9*(a2*a2)*(MW*MW)*(st*st)*A0(Mh1)*(cos(2*theta)*cos(2*theta))*(sqrt(g0sq)*sqrt(g0sq)))/(2.*(g0sq*g0sq)*(Mh2*Mh2)) - 
        (3*a2*(st*st)*A0(Mh2)*(cos(2*theta)*cos(2*theta))*(sqrt(g0sq)*sqrt(g0sq)))/(2.*g0sq) - 
        (3*a2*(Mh1*Mh1)*(st*st)*A0(Mh2)*(cos(2*theta)*cos(2*theta))*(sqrt(g0sq)*sqrt(g0sq)))/(4.*g0sq*(Mh2*Mh2)) + 
        (3*(Mh1*Mh1)*(st*st)*A0(Mh2)*(cos(2*theta)*cos(2*theta))*(sqrt(g0sq)*sqrt(g0sq)))/(16.*(MW*MW)) + 
        (3*(Mh2*Mh2)*(st*st)*A0(Mh2)*(cos(2*theta)*cos(2*theta))*(sqrt(g0sq)*sqrt(g0sq)))/(32.*(MW*MW)) + 
        (9*(a2*a2)*(MW*MW)*(st*st)*A0(Mh2)*(cos(2*theta)*cos(2*theta))*(sqrt(g0sq)*sqrt(g0sq)))/(2.*(g0sq*g0sq)*(Mh2*Mh2)) + 
        (9*a2*ct*A0(Mh1)*cos(3*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(8.*g0sq) - 
        (27*ct*(Mh1*Mh1)*A0(Mh1)*cos(3*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(64.*(MW*MW)) + 
        (9*(a2*a2)*ct*(MW*MW)*A0(Mh1)*cos(3*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(4.*(g0sq*g0sq)*(Mh1*Mh1)) + 
        (3*b3*(st*st*st)*A0(Mh1)*cos(3*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(4.*sqrt(g0sq)*MW) - 
        (3*a2*b3*MW*(st*st*st)*A0(Mh1)*cos(3*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(g0sq*sqrt(g0sq)*(Mh1*Mh1)) + 
        (3*a2*ct*A0(Mh2)*cos(3*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(8.*g0sq) + 
        (3*a2*ct*(Mh2*Mh2)*A0(Mh2)*cos(3*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(8.*g0sq*(Mh1*Mh1)) - 
        (3*ct*(Mh1*Mh1)*A0(Mh2)*cos(3*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(64.*(MW*MW)) - 
        (3*ct*(Mh2*Mh2)*A0(Mh2)*cos(3*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(32.*(MW*MW)) - 
        (3*(a2*a2)*ct*(MW*MW)*A0(Mh2)*cos(3*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(4.*(g0sq*g0sq)*(Mh1*Mh1)) - 
        (b3*st*A0(Mh2)*cos(3*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(16.*sqrt(g0sq)*MW) - 
        (b3*(Mh1*Mh1)*st*A0(Mh2)*cos(3*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(8.*sqrt(g0sq)*(Mh2*Mh2)*MW) + 
        (a2*b3*MW*st*A0(Mh2)*cos(3*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(4.*g0sq*sqrt(g0sq)*(Mh2*Mh2)) - 
        (9*a2*ct*(Mt*Mt)*A0(Mt)*cos(3*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(g0sq*(Mh1*Mh1)) + 
        (9*ct*(Mt*Mt)*A0(Mt)*cos(3*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(4.*(MW*MW)) + 
        (3*a2*ct*A0(MW)*cos(3*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(4.*g0sq) - 
        (3*ct*(Mh1*Mh1)*A0(MW)*cos(3*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(16.*(MW*MW)) + 
        (3*a2*ct*A0(MZ)*cos(3*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(8.*g0sq) - 
        (3*ct*(Mh1*Mh1)*A0(MZ)*cos(3*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(32.*(MW*MW)) - 
        (3*a2*ct*A0(Mh2)*cos(2*theta)*cos(3*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(4.*g0sq) - 
        (3*a2*ct*(Mh2*Mh2)*A0(Mh2)*cos(2*theta)*cos(3*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(8.*g0sq*(Mh1*Mh1)) + 
        (3*ct*(Mh1*Mh1)*A0(Mh2)*cos(2*theta)*cos(3*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(64.*(MW*MW)) + 
        (3*ct*(Mh2*Mh2)*A0(Mh2)*cos(2*theta)*cos(3*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(32.*(MW*MW)) + 
        (9*(a2*a2)*ct*(MW*MW)*A0(Mh2)*cos(2*theta)*cos(3*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(4.*(g0sq*g0sq)*(Mh1*Mh1)) - 
        (b3*st*A0(Mh2)*cos(2*theta)*cos(3*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(16.*sqrt(g0sq)*MW) - 
        (b3*(Mh1*Mh1)*st*A0(Mh2)*cos(2*theta)*cos(3*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(8.*sqrt(g0sq)*(Mh2*Mh2)*MW) + 
        (3*a2*b3*MW*st*A0(Mh2)*cos(2*theta)*cos(3*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(4.*g0sq*sqrt(g0sq)*(Mh2*Mh2)) + 
        (9*a2*A0(Mh1)*(cos(3*theta)*cos(3*theta))*(sqrt(g0sq)*sqrt(g0sq)))/(16.*g0sq) - 
        (9*(Mh1*Mh1)*A0(Mh1)*(cos(3*theta)*cos(3*theta))*(sqrt(g0sq)*sqrt(g0sq)))/(128.*(MW*MW)) - 
        (9*(a2*a2)*(MW*MW)*A0(Mh1)*(cos(3*theta)*cos(3*theta))*(sqrt(g0sq)*sqrt(g0sq)))/(8.*(g0sq*g0sq)*(Mh1*Mh1)) - 
        (b3*(st*st)*A0(Mh1)*sin(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(4.*sqrt(g0sq)*MW) - 
        (b3*(Mh1*Mh1)*(st*st)*A0(Mh1)*sin(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(2.*sqrt(g0sq)*(Mh2*Mh2)*MW) + 
        (a2*b3*MW*(st*st)*A0(Mh1)*sin(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(g0sq*sqrt(g0sq)*(Mh2*Mh2)) + 
        (9*b3*(ct*ct)*A0(Mh2)*sin(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(16.*sqrt(g0sq)*MW) + 
        (3*a2*b3*(ct*ct)*MW*A0(Mh2)*sin(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(4.*g0sq*sqrt(g0sq)*(Mh1*Mh1)) - 
        (3*(b3*b3)*ct*st*A0(Mh2)*sin(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(4.*g0sq*(Mh2*Mh2)) - 
        (3*b3*(st*st)*A0(Mh2)*sin(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(8.*sqrt(g0sq)*MW) - 
        (3*a2*b3*MW*(st*st)*A0(Mh2)*sin(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(2.*g0sq*sqrt(g0sq)*(Mh2*Mh2)) - 
        (b3*b3*ct*(st*st*st)*A0(Mh2)*sin(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(g0sq*(Mh1*Mh1)) + 
        (6*b3*(Mt*Mt)*(st*st)*A0(Mt)*sin(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(sqrt(g0sq)*(Mh2*Mh2)*MW) - 
        (b3*(st*st)*A0(MW)*sin(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(2.*sqrt(g0sq)*MW) - 
        (b3*(st*st)*A0(MZ)*sin(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(4.*sqrt(g0sq)*MW) - 
        (b3*(st*st)*A0(Mh1)*cos(2*theta)*sin(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(4.*sqrt(g0sq)*MW) - 
        (b3*(Mh1*Mh1)*(st*st)*A0(Mh1)*cos(2*theta)*sin(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(2.*sqrt(g0sq)*(Mh2*Mh2)*MW) + 
        (3*a2*b3*MW*(st*st)*A0(Mh1)*cos(2*theta)*sin(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(g0sq*sqrt(g0sq)*(Mh2*Mh2)) + 
        (3*b3*(st*st)*A0(Mh2)*cos(2*theta)*sin(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(8.*sqrt(g0sq)*MW) - 
        (3*a2*b3*MW*(st*st)*A0(Mh2)*cos(2*theta)*sin(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(2.*g0sq*sqrt(g0sq)*(Mh2*Mh2)) + 
        (3*b3*ct*A0(Mh2)*cos(3*theta)*sin(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(16.*sqrt(g0sq)*MW) - 
        (3*a2*b3*ct*MW*A0(Mh2)*cos(3*theta)*sin(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(4.*g0sq*sqrt(g0sq)*(Mh1*Mh1)) - 
        (b3*b3*st*A0(Mh2)*cos(3*theta)*sin(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(4.*g0sq*(Mh2*Mh2)) - 
        (b3*b3*(st*st)*A0(Mh1)*(sin(2*theta)*sin(2*theta))*(sqrt(g0sq)*sqrt(g0sq)))/(2.*g0sq*(Mh2*Mh2))
	);

	Float res = others.fin() + epsPart.div();
    return (double) res;
}

double Renormalization::CalcSelfEnergyH2(double k) {

	// Part proportional to epsilon
	Integral epsPart = (
        -1.0*(g0sq*(st*st)*A0(MW)) - (g0sq*(MZ*MZ)*(st*st)*A0(MZ))/(2.*(MW*MW)) - 2*g0sq*(MW*MW)*(st*st)*B0(k,MW,MW) - 
        (g0sq*(MZ*MZ*MZ*MZ)*(st*st)*B0(k,MZ,MZ))/(MW*MW) + (ct*ct*sqrt(g0sq)*A0(MW)*sqrt(g0sq))/2. + 
        (ct*ct*sqrt(g0sq)*(Mh2*Mh2)*A0(MW)*sqrt(g0sq))/(Mh1*Mh1) - 
        (2*a2*(ct*ct)*(MW*MW)*A0(MW)*sqrt(g0sq))/(sqrt(g0sq)*(Mh1*Mh1)) + (3*b3*ct*MW*st*A0(MW)*sqrt(g0sq))/(Mh2*Mh2) + 
        (3*sqrt(g0sq)*(st*st)*A0(MW)*sqrt(g0sq))/2. + (6*a2*(MW*MW)*(st*st)*A0(MW)*sqrt(g0sq))/(sqrt(g0sq)*(Mh2*Mh2)) - 
        (a2*(ct*ct)*(MZ*MZ)*A0(MZ)*sqrt(g0sq))/(sqrt(g0sq)*(Mh1*Mh1)) + 
        (ct*ct*sqrt(g0sq)*(MZ*MZ)*A0(MZ)*sqrt(g0sq))/(4.*(MW*MW)) + 
        (ct*ct*sqrt(g0sq)*(Mh2*Mh2)*(MZ*MZ)*A0(MZ)*sqrt(g0sq))/(2.*(Mh1*Mh1)*(MW*MW)) + 
        (3*b3*ct*(MZ*MZ)*st*A0(MZ)*sqrt(g0sq))/(2.*(Mh2*Mh2)*MW) + 
        (3*a2*(MZ*MZ)*(st*st)*A0(MZ)*sqrt(g0sq))/(sqrt(g0sq)*(Mh2*Mh2)) + 
        (3*sqrt(g0sq)*(MZ*MZ)*(st*st)*A0(MZ)*sqrt(g0sq))/(4.*(MW*MW)) - 
        (ct*ct*sqrt(g0sq)*A0(MW)*cos(2*theta)*sqrt(g0sq))/2. - 
        (ct*ct*sqrt(g0sq)*(Mh2*Mh2)*A0(MW)*cos(2*theta)*sqrt(g0sq))/(Mh1*Mh1) + 
        (6*a2*(ct*ct)*(MW*MW)*A0(MW)*cos(2*theta)*sqrt(g0sq))/(sqrt(g0sq)*(Mh1*Mh1)) - 
        (3*sqrt(g0sq)*(st*st)*A0(MW)*cos(2*theta)*sqrt(g0sq))/2. + 
        (6*a2*(MW*MW)*(st*st)*A0(MW)*cos(2*theta)*sqrt(g0sq))/(sqrt(g0sq)*(Mh2*Mh2)) + 
        (3*a2*(ct*ct)*(MZ*MZ)*A0(MZ)*cos(2*theta)*sqrt(g0sq))/(sqrt(g0sq)*(Mh1*Mh1)) - 
        (ct*ct*sqrt(g0sq)*(MZ*MZ)*A0(MZ)*cos(2*theta)*sqrt(g0sq))/(4.*(MW*MW)) - 
        (ct*ct*sqrt(g0sq)*(Mh2*Mh2)*(MZ*MZ)*A0(MZ)*cos(2*theta)*sqrt(g0sq))/(2.*(Mh1*Mh1)*(MW*MW)) + 
        (3*a2*(MZ*MZ)*(st*st)*A0(MZ)*cos(2*theta)*sqrt(g0sq))/(sqrt(g0sq)*(Mh2*Mh2)) - 
        (3*sqrt(g0sq)*(MZ*MZ)*(st*st)*A0(MZ)*cos(2*theta)*sqrt(g0sq))/(4.*(MW*MW)) + 
        (b3*MW*st*A0(MW)*cos(3*theta)*sqrt(g0sq))/(Mh2*Mh2) + 
        (b3*(MZ*MZ)*st*A0(MZ)*cos(3*theta)*sqrt(g0sq))/(2.*(Mh2*Mh2)*MW) - 
        (2*b3*(ct*ct)*MW*A0(MW)*sin(2*theta)*sqrt(g0sq))/(Mh1*Mh1) - 
        (b3*(ct*ct)*(MZ*MZ)*A0(MZ)*sin(2*theta)*sqrt(g0sq))/(Mh1*Mh1*MW)
	);

	Integral others = (
        (a2*A0(Mh1))/8. + (3*b4*A0(Mh1))/8. + (3*g0sq*(Mh1*Mh1)*A0(Mh1))/(128.*(MW*MW)) + 
        (3*g0sq*(Mh2*Mh2)*A0(Mh1))/(128.*(MW*MW)) + 3*b4*(ct*ct*ct*ct)*A0(Mh2) + 3*a2*(ct*ct)*(st*st)*A0(Mh2) + 
        (3*g0sq*(Mh1*Mh1)*(st*st*st*st)*A0(Mh2))/(16.*(MW*MW)) + (3*g0sq*(Mh2*Mh2)*(st*st*st*st)*A0(Mh2))/(16.*(MW*MW)) + 
        a2*(ct*ct)*A0(MW) + 2*g0sq*(st*st)*A0(MW) + (g0sq*(Mh1*Mh1)*(st*st)*A0(MW))/(8.*(MW*MW)) + 
        (g0sq*(Mh2*Mh2)*(st*st)*A0(MW))/(8.*(MW*MW)) + (a2*(ct*ct)*A0(MZ))/2. + 
        (g0sq*(Mh1*Mh1)*(st*st)*A0(MZ))/(16.*(MW*MW)) + (g0sq*(Mh2*Mh2)*(st*st)*A0(MZ))/(16.*(MW*MW)) + 
        (g0sq*(MZ*MZ)*(st*st)*A0(MZ))/(MW*MW) - (a2*(Mh1*Mh1)*(st*st)*B0(k,Mh1,Mh1))/2. - 
        (a2*(Mh2*Mh2)*(st*st)*B0(k,Mh1,Mh1))/4. + (g0sq*(Mh1*Mh1*Mh1*Mh1)*(st*st)*B0(k,Mh1,Mh1))/(8.*(MW*MW)) + 
        (g0sq*(Mh1*Mh1)*(Mh2*Mh2)*(st*st)*B0(k,Mh1,Mh1))/(8.*(MW*MW)) + 
        (g0sq*(Mh2*Mh2*Mh2*Mh2)*(st*st)*B0(k,Mh1,Mh1))/(32.*(MW*MW)) + (a2*a2*(MW*MW)*(st*st)*B0(k,Mh1,Mh1))/(2.*g0sq) - 
        (a2*(ct*ct)*(Mh1*Mh1)*B0(k,Mh1,Mh2))/2. - a2*(ct*ct)*(Mh2*Mh2)*B0(k,Mh1,Mh2) + 
        (ct*ct*g0sq*(Mh1*Mh1*Mh1*Mh1)*B0(k,Mh1,Mh2))/(16.*(MW*MW)) + 
        (ct*ct*g0sq*(Mh1*Mh1)*(Mh2*Mh2)*B0(k,Mh1,Mh2))/(4.*(MW*MW)) + 
        (ct*ct*g0sq*(Mh2*Mh2*Mh2*Mh2)*B0(k,Mh1,Mh2))/(4.*(MW*MW)) + (a2*a2*(ct*ct)*(MW*MW)*B0(k,Mh1,Mh2))/g0sq + 
        (9*(b3*b3)*(ct*ct)*B0(k,Mh2,Mh2))/8. + (9*b3*ct*sqrt(g0sq)*(Mh2*Mh2)*st*B0(k,Mh2,Mh2))/(8.*MW) + 
        (9*a2*b3*ct*MW*st*B0(k,Mh2,Mh2))/(2.*sqrt(g0sq)) + (9*a2*(Mh2*Mh2)*(st*st)*B0(k,Mh2,Mh2))/4. + 
        (9*g0sq*(Mh2*Mh2*Mh2*Mh2)*(st*st)*B0(k,Mh2,Mh2))/(32.*(MW*MW)) + 
        (9*(a2*a2)*(MW*MW)*(st*st)*B0(k,Mh2,Mh2))/(2.*g0sq) + (g0sq*(Mh2*Mh2*Mh2*Mh2)*(st*st)*B0(k,MW,MW))/(4.*(MW*MW)) + 
        (7*g0sq*(MW*MW)*(st*st)*B0(k,MW,MW))/2. + (g0sq*(Mh2*Mh2*Mh2*Mh2)*(st*st)*B0(k,MZ,MZ))/(8.*(MW*MW)) + 
        (7*g0sq*(MZ*MZ*MZ*MZ)*(st*st)*B0(k,MZ,MZ))/(4.*(MW*MW)) + (3*g0sq*(Mh1*Mh1)*A0(Mh1)*cos(2*theta))/(256.*(MW*MW)) - 
        (3*g0sq*(Mh2*Mh2)*A0(Mh1)*cos(2*theta))/(256.*(MW*MW)) + 
        (3*g0sq*(Mh1*Mh1)*(st*st*st*st)*A0(Mh2)*cos(2*theta))/(16.*(MW*MW)) - 
        (3*g0sq*(Mh2*Mh2)*(st*st*st*st)*A0(Mh2)*cos(2*theta))/(16.*(MW*MW)) + 
        (g0sq*(Mh1*Mh1)*(st*st)*A0(MW)*cos(2*theta))/(8.*(MW*MW)) - 
        (g0sq*(Mh2*Mh2)*(st*st)*A0(MW)*cos(2*theta))/(8.*(MW*MW)) + 
        (g0sq*(Mh1*Mh1)*(st*st)*A0(MZ)*cos(2*theta))/(16.*(MW*MW)) - 
        (g0sq*(Mh2*Mh2)*(st*st)*A0(MZ)*cos(2*theta))/(16.*(MW*MW)) - 2*a2*(Mh1*Mh1)*(st*st)*B0(k,Mh1,Mh1)*cos(2*theta) - 
        a2*(Mh2*Mh2)*(st*st)*B0(k,Mh1,Mh1)*cos(2*theta) + 
        (g0sq*(Mh1*Mh1*Mh1*Mh1)*(st*st)*B0(k,Mh1,Mh1)*cos(2*theta))/(4.*(MW*MW)) + 
        (g0sq*(Mh1*Mh1)*(Mh2*Mh2)*(st*st)*B0(k,Mh1,Mh1)*cos(2*theta))/(4.*(MW*MW)) + 
        (g0sq*(Mh2*Mh2*Mh2*Mh2)*(st*st)*B0(k,Mh1,Mh1)*cos(2*theta))/(16.*(MW*MW)) + 
        (3*(a2*a2)*(MW*MW)*(st*st)*B0(k,Mh1,Mh1)*cos(2*theta))/g0sq + 2*a2*(ct*ct)*(Mh1*Mh1)*B0(k,Mh1,Mh2)*cos(2*theta) + 
        4*a2*(ct*ct)*(Mh2*Mh2)*B0(k,Mh1,Mh2)*cos(2*theta) - 
        (ct*ct*g0sq*(Mh1*Mh1*Mh1*Mh1)*B0(k,Mh1,Mh2)*cos(2*theta))/(8.*(MW*MW)) - 
        (ct*ct*g0sq*(Mh1*Mh1)*(Mh2*Mh2)*B0(k,Mh1,Mh2)*cos(2*theta))/(2.*(MW*MW)) - 
        (ct*ct*g0sq*(Mh2*Mh2*Mh2*Mh2)*B0(k,Mh1,Mh2)*cos(2*theta))/(2.*(MW*MW)) - 
        (6*(a2*a2)*(ct*ct)*(MW*MW)*B0(k,Mh1,Mh2)*cos(2*theta))/g0sq - 
        (9*b3*ct*sqrt(g0sq)*(Mh2*Mh2)*st*B0(k,Mh2,Mh2)*cos(2*theta))/(8.*MW) + 
        (9*a2*b3*ct*MW*st*B0(k,Mh2,Mh2)*cos(2*theta))/(2.*sqrt(g0sq)) - 
        (9*g0sq*(Mh2*Mh2*Mh2*Mh2)*(st*st)*B0(k,Mh2,Mh2)*cos(2*theta))/(16.*(MW*MW)) + 
        (9*(a2*a2)*(MW*MW)*(st*st)*B0(k,Mh2,Mh2)*cos(2*theta))/g0sq - 
        (3*a2*(Mh1*Mh1)*(st*st)*B0(k,Mh1,Mh1)*(cos(2*theta)*cos(2*theta)))/2. - 
        (3*a2*(Mh2*Mh2)*(st*st)*B0(k,Mh1,Mh1)*(cos(2*theta)*cos(2*theta)))/4. + 
        (g0sq*(Mh1*Mh1*Mh1*Mh1)*(st*st)*B0(k,Mh1,Mh1)*(cos(2*theta)*cos(2*theta)))/(8.*(MW*MW)) + 
        (g0sq*(Mh1*Mh1)*(Mh2*Mh2)*(st*st)*B0(k,Mh1,Mh1)*(cos(2*theta)*cos(2*theta)))/(8.*(MW*MW)) + 
        (g0sq*(Mh2*Mh2*Mh2*Mh2)*(st*st)*B0(k,Mh1,Mh1)*(cos(2*theta)*cos(2*theta)))/(32.*(MW*MW)) + 
        (9*(a2*a2)*(MW*MW)*(st*st)*B0(k,Mh1,Mh1)*(cos(2*theta)*cos(2*theta)))/(2.*g0sq) - 
        (3*a2*(ct*ct)*(Mh1*Mh1)*B0(k,Mh1,Mh2)*(cos(2*theta)*cos(2*theta)))/2. - 
        3*a2*(ct*ct)*(Mh2*Mh2)*B0(k,Mh1,Mh2)*(cos(2*theta)*cos(2*theta)) + 
        (ct*ct*g0sq*(Mh1*Mh1*Mh1*Mh1)*B0(k,Mh1,Mh2)*(cos(2*theta)*cos(2*theta)))/(16.*(MW*MW)) + 
        (ct*ct*g0sq*(Mh1*Mh1)*(Mh2*Mh2)*B0(k,Mh1,Mh2)*(cos(2*theta)*cos(2*theta)))/(4.*(MW*MW)) + 
        (ct*ct*g0sq*(Mh2*Mh2*Mh2*Mh2)*B0(k,Mh1,Mh2)*(cos(2*theta)*cos(2*theta)))/(4.*(MW*MW)) + 
        (9*(a2*a2)*(ct*ct)*(MW*MW)*B0(k,Mh1,Mh2)*(cos(2*theta)*cos(2*theta)))/g0sq - 
        (9*a2*(Mh2*Mh2)*(st*st)*B0(k,Mh2,Mh2)*(cos(2*theta)*cos(2*theta)))/4. + 
        (9*g0sq*(Mh2*Mh2*Mh2*Mh2)*(st*st)*B0(k,Mh2,Mh2)*(cos(2*theta)*cos(2*theta)))/(32.*(MW*MW)) + 
        (9*(a2*a2)*(MW*MW)*(st*st)*B0(k,Mh2,Mh2)*(cos(2*theta)*cos(2*theta)))/(2.*g0sq) + 
        (3*(b3*b3)*ct*B0(k,Mh2,Mh2)*cos(3*theta))/4. + (3*b3*sqrt(g0sq)*(Mh2*Mh2)*st*B0(k,Mh2,Mh2)*cos(3*theta))/(8.*MW) + 
        (3*a2*b3*MW*st*B0(k,Mh2,Mh2)*cos(3*theta))/(2.*sqrt(g0sq)) - 
        (3*b3*sqrt(g0sq)*(Mh2*Mh2)*st*B0(k,Mh2,Mh2)*cos(2*theta)*cos(3*theta))/(8.*MW) + 
        (3*a2*b3*MW*st*B0(k,Mh2,Mh2)*cos(2*theta)*cos(3*theta))/(2.*sqrt(g0sq)) + 
        (b3*b3*B0(k,Mh2,Mh2)*(cos(3*theta)*cos(3*theta)))/8. + (3*a2*A0(Mh1)*cos(4*theta))/8. - 
        (3*b4*A0(Mh1)*cos(4*theta))/8. - (3*g0sq*(Mh1*Mh1)*A0(Mh1)*cos(4*theta))/(128.*(MW*MW)) - 
        (3*g0sq*(Mh2*Mh2)*A0(Mh1)*cos(4*theta))/(128.*(MW*MW)) - (3*g0sq*(Mh1*Mh1)*A0(Mh1)*cos(6*theta))/(256.*(MW*MW)) + 
        (3*g0sq*(Mh2*Mh2)*A0(Mh1)*cos(6*theta))/(256.*(MW*MW)) - (3*g0sq*(Mt*Mt)*(st*st)*K(k,Mt,Mt,1,1,Mt*Mt))/(MW*MW) - 
        (g0sq*(st*st)*K(k,MW,MW,1,4,4*(k*k)))/2. - (g0sq*(MZ*MZ)*(st*st)*K(k,MZ,MZ,1,4,4*(k*k)))/(4.*(MW*MW)) + 
        (b3*sqrt(g0sq)*(Mh1*Mh1)*(st*st)*B0(k,Mh1,Mh1)*sin(2*theta))/(2.*MW) + 
        (b3*sqrt(g0sq)*(Mh2*Mh2)*(st*st)*B0(k,Mh1,Mh1)*sin(2*theta))/(4.*MW) - 
        (a2*b3*MW*(st*st)*B0(k,Mh1,Mh1)*sin(2*theta))/sqrt(g0sq) - 
        (b3*(ct*ct)*sqrt(g0sq)*(Mh1*Mh1)*B0(k,Mh1,Mh2)*sin(2*theta))/(2.*MW) - 
        (b3*(ct*ct)*sqrt(g0sq)*(Mh2*Mh2)*B0(k,Mh1,Mh2)*sin(2*theta))/MW + 
        (2*a2*b3*(ct*ct)*MW*B0(k,Mh1,Mh2)*sin(2*theta))/sqrt(g0sq) + 
        (b3*sqrt(g0sq)*(Mh1*Mh1)*(st*st)*B0(k,Mh1,Mh1)*cos(2*theta)*sin(2*theta))/(2.*MW) + 
        (b3*sqrt(g0sq)*(Mh2*Mh2)*(st*st)*B0(k,Mh1,Mh1)*cos(2*theta)*sin(2*theta))/(4.*MW) - 
        (3*a2*b3*MW*(st*st)*B0(k,Mh1,Mh1)*cos(2*theta)*sin(2*theta))/sqrt(g0sq) + 
        (b3*(ct*ct)*sqrt(g0sq)*(Mh1*Mh1)*B0(k,Mh1,Mh2)*cos(2*theta)*sin(2*theta))/(2.*MW) + 
        (b3*(ct*ct)*sqrt(g0sq)*(Mh2*Mh2)*B0(k,Mh1,Mh2)*cos(2*theta)*sin(2*theta))/MW - 
        (6*a2*b3*(ct*ct)*MW*B0(k,Mh1,Mh2)*cos(2*theta)*sin(2*theta))/sqrt(g0sq) + 
        (b3*b3*(st*st)*B0(k,Mh1,Mh1)*(sin(2*theta)*sin(2*theta)))/2. + 
        b3*b3*(ct*ct)*B0(k,Mh1,Mh2)*(sin(2*theta)*sin(2*theta)) - (3*(ct*ct)*sqrt(g0sq)*A0(MW)*sqrt(g0sq))/4. - 
        (3*(ct*ct)*sqrt(g0sq)*(Mh2*Mh2)*A0(MW)*sqrt(g0sq))/(2.*(Mh1*Mh1)) + 
        (3*a2*(ct*ct)*(MW*MW)*A0(MW)*sqrt(g0sq))/(sqrt(g0sq)*(Mh1*Mh1)) - (9*b3*ct*MW*st*A0(MW)*sqrt(g0sq))/(2.*(Mh2*Mh2)) - 
        (9*sqrt(g0sq)*(st*st)*A0(MW)*sqrt(g0sq))/4. - (9*a2*(MW*MW)*(st*st)*A0(MW)*sqrt(g0sq))/(sqrt(g0sq)*(Mh2*Mh2)) + 
        (3*a2*(ct*ct)*(MZ*MZ)*A0(MZ)*sqrt(g0sq))/(2.*sqrt(g0sq)*(Mh1*Mh1)) - 
        (3*(ct*ct)*sqrt(g0sq)*(MZ*MZ)*A0(MZ)*sqrt(g0sq))/(8.*(MW*MW)) - 
        (3*(ct*ct)*sqrt(g0sq)*(Mh2*Mh2)*(MZ*MZ)*A0(MZ)*sqrt(g0sq))/(4.*(Mh1*Mh1)*(MW*MW)) - 
        (9*b3*ct*(MZ*MZ)*st*A0(MZ)*sqrt(g0sq))/(4.*(Mh2*Mh2)*MW) - 
        (9*a2*(MZ*MZ)*(st*st)*A0(MZ)*sqrt(g0sq))/(2.*sqrt(g0sq)*(Mh2*Mh2)) - 
        (9*sqrt(g0sq)*(MZ*MZ)*(st*st)*A0(MZ)*sqrt(g0sq))/(8.*(MW*MW)) + 
        (3*(ct*ct)*sqrt(g0sq)*A0(MW)*cos(2*theta)*sqrt(g0sq))/4. + 
        (3*(ct*ct)*sqrt(g0sq)*(Mh2*Mh2)*A0(MW)*cos(2*theta)*sqrt(g0sq))/(2.*(Mh1*Mh1)) - 
        (9*a2*(ct*ct)*(MW*MW)*A0(MW)*cos(2*theta)*sqrt(g0sq))/(sqrt(g0sq)*(Mh1*Mh1)) + 
        (9*sqrt(g0sq)*(st*st)*A0(MW)*cos(2*theta)*sqrt(g0sq))/4. - 
        (9*a2*(MW*MW)*(st*st)*A0(MW)*cos(2*theta)*sqrt(g0sq))/(sqrt(g0sq)*(Mh2*Mh2)) - 
        (9*a2*(ct*ct)*(MZ*MZ)*A0(MZ)*cos(2*theta)*sqrt(g0sq))/(2.*sqrt(g0sq)*(Mh1*Mh1)) + 
        (3*(ct*ct)*sqrt(g0sq)*(MZ*MZ)*A0(MZ)*cos(2*theta)*sqrt(g0sq))/(8.*(MW*MW)) + 
        (3*(ct*ct)*sqrt(g0sq)*(Mh2*Mh2)*(MZ*MZ)*A0(MZ)*cos(2*theta)*sqrt(g0sq))/(4.*(Mh1*Mh1)*(MW*MW)) - 
        (9*a2*(MZ*MZ)*(st*st)*A0(MZ)*cos(2*theta)*sqrt(g0sq))/(2.*sqrt(g0sq)*(Mh2*Mh2)) + 
        (9*sqrt(g0sq)*(MZ*MZ)*(st*st)*A0(MZ)*cos(2*theta)*sqrt(g0sq))/(8.*(MW*MW)) - 
        (3*b3*MW*st*A0(MW)*cos(3*theta)*sqrt(g0sq))/(2.*(Mh2*Mh2)) - 
        (3*b3*(MZ*MZ)*st*A0(MZ)*cos(3*theta)*sqrt(g0sq))/(4.*(Mh2*Mh2)*MW) + 
        (3*b3*(ct*ct)*MW*A0(MW)*sin(2*theta)*sqrt(g0sq))/(Mh1*Mh1) + 
        (3*b3*(ct*ct)*(MZ*MZ)*A0(MZ)*sin(2*theta)*sqrt(g0sq))/(2.*(Mh1*Mh1)*MW) + 
        (3*a2*(ct*ct)*A0(Mh1)*(sqrt(g0sq)*sqrt(g0sq)))/(8.*g0sq) - 
        (3*a2*(ct*ct)*(Mh2*Mh2)*A0(Mh1)*(sqrt(g0sq)*sqrt(g0sq)))/(8.*g0sq*(Mh1*Mh1)) - 
        (9*(ct*ct)*(Mh1*Mh1)*A0(Mh1)*(sqrt(g0sq)*sqrt(g0sq)))/(64.*(MW*MW)) - 
        (9*(ct*ct)*(Mh2*Mh2)*A0(Mh1)*(sqrt(g0sq)*sqrt(g0sq)))/(32.*(MW*MW)) + 
        (3*(a2*a2)*(ct*ct)*(MW*MW)*A0(Mh1)*(sqrt(g0sq)*sqrt(g0sq)))/(4.*(g0sq*g0sq)*(Mh1*Mh1)) - 
        (3*b3*ct*st*A0(Mh1)*(sqrt(g0sq)*sqrt(g0sq)))/(16.*sqrt(g0sq)*MW) - 
        (3*b3*ct*(Mh1*Mh1)*st*A0(Mh1)*(sqrt(g0sq)*sqrt(g0sq)))/(8.*sqrt(g0sq)*(Mh2*Mh2)*MW) + 
        (3*a2*b3*ct*MW*st*A0(Mh1)*(sqrt(g0sq)*sqrt(g0sq)))/(4.*g0sq*sqrt(g0sq)*(Mh2*Mh2)) - 
        (3*a2*(Mh1*Mh1)*(st*st)*A0(Mh1)*(sqrt(g0sq)*sqrt(g0sq)))/(4.*g0sq*(Mh2*Mh2)) - 
        (3*(Mh1*Mh1)*(st*st)*A0(Mh1)*(sqrt(g0sq)*sqrt(g0sq)))/(16.*(MW*MW)) - 
        (3*(Mh2*Mh2)*(st*st)*A0(Mh1)*(sqrt(g0sq)*sqrt(g0sq)))/(32.*(MW*MW)) + 
        (3*(a2*a2)*(MW*MW)*(st*st)*A0(Mh1)*(sqrt(g0sq)*sqrt(g0sq)))/(2.*(g0sq*g0sq)*(Mh2*Mh2)) + 
        (b3*ct*(st*st*st)*A0(Mh1)*(sqrt(g0sq)*sqrt(g0sq)))/(4.*sqrt(g0sq)*MW) + 
        (b3*ct*(Mh2*Mh2)*(st*st*st)*A0(Mh1)*(sqrt(g0sq)*sqrt(g0sq)))/(2.*sqrt(g0sq)*(Mh1*Mh1)*MW) - 
        (a2*b3*ct*MW*(st*st*st)*A0(Mh1)*(sqrt(g0sq)*sqrt(g0sq)))/(g0sq*sqrt(g0sq)*(Mh1*Mh1)) + 
        (a2*(ct*ct)*A0(Mh2)*(sqrt(g0sq)*sqrt(g0sq)))/(4.*g0sq) - 
        (9*(b3*b3)*(ct*ct)*A0(Mh2)*(sqrt(g0sq)*sqrt(g0sq)))/(8.*g0sq*(Mh2*Mh2)) + 
        (a2*(ct*ct)*(Mh2*Mh2)*A0(Mh2)*(sqrt(g0sq)*sqrt(g0sq)))/(2.*g0sq*(Mh1*Mh1)) - 
        (ct*ct*(Mh1*Mh1)*A0(Mh2)*(sqrt(g0sq)*sqrt(g0sq)))/(32.*(MW*MW)) - 
        (ct*ct*(Mh2*Mh2)*A0(Mh2)*(sqrt(g0sq)*sqrt(g0sq)))/(8.*(MW*MW)) - 
        (ct*ct*(Mh2*Mh2*Mh2*Mh2)*A0(Mh2)*(sqrt(g0sq)*sqrt(g0sq)))/(8.*(Mh1*Mh1)*(MW*MW)) - 
        (a2*a2*(ct*ct)*(MW*MW)*A0(Mh2)*(sqrt(g0sq)*sqrt(g0sq)))/(2.*(g0sq*g0sq)*(Mh1*Mh1)) - 
        (9*b3*ct*st*A0(Mh2)*(sqrt(g0sq)*sqrt(g0sq)))/(8.*sqrt(g0sq)*MW) - 
        (9*a2*b3*ct*MW*st*A0(Mh2)*(sqrt(g0sq)*sqrt(g0sq)))/(2.*g0sq*sqrt(g0sq)*(Mh2*Mh2)) - 
        (9*a2*(st*st)*A0(Mh2)*(sqrt(g0sq)*sqrt(g0sq)))/(4.*g0sq) - 
        (9*(Mh2*Mh2)*(st*st)*A0(Mh2)*(sqrt(g0sq)*sqrt(g0sq)))/(32.*(MW*MW)) - 
        (9*(a2*a2)*(MW*MW)*(st*st)*A0(Mh2)*(sqrt(g0sq)*sqrt(g0sq)))/(2.*(g0sq*g0sq)*(Mh2*Mh2)) - 
        (6*a2*(ct*ct)*(Mt*Mt)*A0(Mt)*(sqrt(g0sq)*sqrt(g0sq)))/(g0sq*(Mh1*Mh1)) + 
        (3*(ct*ct)*(Mt*Mt)*A0(Mt)*(sqrt(g0sq)*sqrt(g0sq)))/(2.*(MW*MW)) + 
        (3*(ct*ct)*(Mh2*Mh2)*(Mt*Mt)*A0(Mt)*(sqrt(g0sq)*sqrt(g0sq)))/(Mh1*Mh1*(MW*MW)) + 
        (9*b3*ct*(Mt*Mt)*st*A0(Mt)*(sqrt(g0sq)*sqrt(g0sq)))/(sqrt(g0sq)*(Mh2*Mh2)*MW) + 
        (18*a2*(Mt*Mt)*(st*st)*A0(Mt)*(sqrt(g0sq)*sqrt(g0sq)))/(g0sq*(Mh2*Mh2)) + 
        (9*(Mt*Mt)*(st*st)*A0(Mt)*(sqrt(g0sq)*sqrt(g0sq)))/(2.*(MW*MW)) + 
        (a2*(ct*ct)*A0(MW)*(sqrt(g0sq)*sqrt(g0sq)))/(2.*g0sq) - 
        (ct*ct*(Mh1*Mh1)*A0(MW)*(sqrt(g0sq)*sqrt(g0sq)))/(8.*(MW*MW)) - 
        (ct*ct*(Mh2*Mh2)*A0(MW)*(sqrt(g0sq)*sqrt(g0sq)))/(4.*(MW*MW)) - 
        (3*b3*ct*st*A0(MW)*(sqrt(g0sq)*sqrt(g0sq)))/(4.*sqrt(g0sq)*MW) - 
        (3*a2*(st*st)*A0(MW)*(sqrt(g0sq)*sqrt(g0sq)))/(2.*g0sq) - 
        (3*(Mh2*Mh2)*(st*st)*A0(MW)*(sqrt(g0sq)*sqrt(g0sq)))/(8.*(MW*MW)) + 
        (a2*(ct*ct)*A0(MZ)*(sqrt(g0sq)*sqrt(g0sq)))/(4.*g0sq) - 
        (ct*ct*(Mh1*Mh1)*A0(MZ)*(sqrt(g0sq)*sqrt(g0sq)))/(16.*(MW*MW)) - 
        (ct*ct*(Mh2*Mh2)*A0(MZ)*(sqrt(g0sq)*sqrt(g0sq)))/(8.*(MW*MW)) - 
        (3*b3*ct*st*A0(MZ)*(sqrt(g0sq)*sqrt(g0sq)))/(8.*sqrt(g0sq)*MW) - 
        (3*a2*(st*st)*A0(MZ)*(sqrt(g0sq)*sqrt(g0sq)))/(4.*g0sq) - 
        (3*(Mh2*Mh2)*(st*st)*A0(MZ)*(sqrt(g0sq)*sqrt(g0sq)))/(16.*(MW*MW)) - 
        (3*a2*(ct*ct)*A0(Mh1)*cos(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(2.*g0sq) + 
        (3*a2*(ct*ct)*(Mh2*Mh2)*A0(Mh1)*cos(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(8.*g0sq*(Mh1*Mh1)) + 
        (9*(ct*ct)*(Mh1*Mh1)*A0(Mh1)*cos(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(64.*(MW*MW)) + 
        (9*(ct*ct)*(Mh2*Mh2)*A0(Mh1)*cos(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(32.*(MW*MW)) - 
        (9*(a2*a2)*(ct*ct)*(MW*MW)*A0(Mh1)*cos(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(4.*(g0sq*g0sq)*(Mh1*Mh1)) - 
        (3*b3*ct*st*A0(Mh1)*cos(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(16.*sqrt(g0sq)*MW) - 
        (3*b3*ct*(Mh1*Mh1)*st*A0(Mh1)*cos(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(8.*sqrt(g0sq)*(Mh2*Mh2)*MW) + 
        (9*a2*b3*ct*MW*st*A0(Mh1)*cos(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(4.*g0sq*sqrt(g0sq)*(Mh2*Mh2)) - 
        (3*a2*(Mh1*Mh1)*(st*st)*A0(Mh1)*cos(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(2.*g0sq*(Mh2*Mh2)) + 
        (6*(a2*a2)*(MW*MW)*(st*st)*A0(Mh1)*cos(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(g0sq*g0sq*(Mh2*Mh2)) - 
        (b3*ct*(st*st*st)*A0(Mh1)*cos(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(4.*sqrt(g0sq)*MW) - 
        (b3*ct*(Mh2*Mh2)*(st*st*st)*A0(Mh1)*cos(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(2.*sqrt(g0sq)*(Mh1*Mh1)*MW) + 
        (3*a2*b3*ct*MW*(st*st*st)*A0(Mh1)*cos(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(g0sq*sqrt(g0sq)*(Mh1*Mh1)) - 
        (a2*(ct*ct)*A0(Mh2)*cos(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/g0sq - 
        (2*a2*(ct*ct)*(Mh2*Mh2)*A0(Mh2)*cos(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(g0sq*(Mh1*Mh1)) + 
        (ct*ct*(Mh1*Mh1)*A0(Mh2)*cos(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(16.*(MW*MW)) + 
        (ct*ct*(Mh2*Mh2)*A0(Mh2)*cos(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(4.*(MW*MW)) + 
        (ct*ct*(Mh2*Mh2*Mh2*Mh2)*A0(Mh2)*cos(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(4.*(Mh1*Mh1)*(MW*MW)) + 
        (3*(a2*a2)*(ct*ct)*(MW*MW)*A0(Mh2)*cos(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(g0sq*g0sq*(Mh1*Mh1)) + 
        (9*b3*ct*st*A0(Mh2)*cos(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(8.*sqrt(g0sq)*MW) - 
        (9*a2*b3*ct*MW*st*A0(Mh2)*cos(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(2.*g0sq*sqrt(g0sq)*(Mh2*Mh2)) + 
        (9*(Mh2*Mh2)*(st*st)*A0(Mh2)*cos(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(16.*(MW*MW)) - 
        (9*(a2*a2)*(MW*MW)*(st*st)*A0(Mh2)*cos(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(g0sq*g0sq*(Mh2*Mh2)) + 
        (18*a2*(ct*ct)*(Mt*Mt)*A0(Mt)*cos(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(g0sq*(Mh1*Mh1)) - 
        (3*(ct*ct)*(Mt*Mt)*A0(Mt)*cos(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(2.*(MW*MW)) - 
        (3*(ct*ct)*(Mh2*Mh2)*(Mt*Mt)*A0(Mt)*cos(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(Mh1*Mh1*(MW*MW)) + 
        (18*a2*(Mt*Mt)*(st*st)*A0(Mt)*cos(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(g0sq*(Mh2*Mh2)) - 
        (9*(Mt*Mt)*(st*st)*A0(Mt)*cos(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(2.*(MW*MW)) - 
        (3*a2*(ct*ct)*A0(MW)*cos(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(2.*g0sq) + 
        (ct*ct*(Mh1*Mh1)*A0(MW)*cos(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(8.*(MW*MW)) + 
        (ct*ct*(Mh2*Mh2)*A0(MW)*cos(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(4.*(MW*MW)) - 
        (3*a2*(st*st)*A0(MW)*cos(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(2.*g0sq) + 
        (3*(Mh2*Mh2)*(st*st)*A0(MW)*cos(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(8.*(MW*MW)) - 
        (3*a2*(ct*ct)*A0(MZ)*cos(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(4.*g0sq) + 
        (ct*ct*(Mh1*Mh1)*A0(MZ)*cos(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(16.*(MW*MW)) + 
        (ct*ct*(Mh2*Mh2)*A0(MZ)*cos(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(8.*(MW*MW)) - 
        (3*a2*(st*st)*A0(MZ)*cos(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(4.*g0sq) + 
        (3*(Mh2*Mh2)*(st*st)*A0(MZ)*cos(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(16.*(MW*MW)) - 
        (3*a2*(st*st)*A0(Mh1)*(cos(2*theta)*cos(2*theta))*(sqrt(g0sq)*sqrt(g0sq)))/(2.*g0sq) - 
        (3*a2*(Mh1*Mh1)*(st*st)*A0(Mh1)*(cos(2*theta)*cos(2*theta))*(sqrt(g0sq)*sqrt(g0sq)))/(4.*g0sq*(Mh2*Mh2)) + 
        (3*(Mh1*Mh1)*(st*st)*A0(Mh1)*(cos(2*theta)*cos(2*theta))*(sqrt(g0sq)*sqrt(g0sq)))/(16.*(MW*MW)) + 
        (3*(Mh2*Mh2)*(st*st)*A0(Mh1)*(cos(2*theta)*cos(2*theta))*(sqrt(g0sq)*sqrt(g0sq)))/(32.*(MW*MW)) + 
        (9*(a2*a2)*(MW*MW)*(st*st)*A0(Mh1)*(cos(2*theta)*cos(2*theta))*(sqrt(g0sq)*sqrt(g0sq)))/(2.*(g0sq*g0sq)*(Mh2*Mh2)) + 
        (3*a2*(ct*ct)*A0(Mh2)*(cos(2*theta)*cos(2*theta))*(sqrt(g0sq)*sqrt(g0sq)))/(4.*g0sq) + 
        (3*a2*(ct*ct)*(Mh2*Mh2)*A0(Mh2)*(cos(2*theta)*cos(2*theta))*(sqrt(g0sq)*sqrt(g0sq)))/(2.*g0sq*(Mh1*Mh1)) - 
        (ct*ct*(Mh1*Mh1)*A0(Mh2)*(cos(2*theta)*cos(2*theta))*(sqrt(g0sq)*sqrt(g0sq)))/(32.*(MW*MW)) - 
        (ct*ct*(Mh2*Mh2)*A0(Mh2)*(cos(2*theta)*cos(2*theta))*(sqrt(g0sq)*sqrt(g0sq)))/(8.*(MW*MW)) - 
        (ct*ct*(Mh2*Mh2*Mh2*Mh2)*A0(Mh2)*(cos(2*theta)*cos(2*theta))*(sqrt(g0sq)*sqrt(g0sq)))/(8.*(Mh1*Mh1)*(MW*MW)) - 
        (9*(a2*a2)*(ct*ct)*(MW*MW)*A0(Mh2)*(cos(2*theta)*cos(2*theta))*(sqrt(g0sq)*sqrt(g0sq)))/(2.*(g0sq*g0sq)*(Mh1*Mh1)) + 
        (9*a2*(st*st)*A0(Mh2)*(cos(2*theta)*cos(2*theta))*(sqrt(g0sq)*sqrt(g0sq)))/(4.*g0sq) - 
        (9*(Mh2*Mh2)*(st*st)*A0(Mh2)*(cos(2*theta)*cos(2*theta))*(sqrt(g0sq)*sqrt(g0sq)))/(32.*(MW*MW)) - 
        (9*(a2*a2)*(MW*MW)*(st*st)*A0(Mh2)*(cos(2*theta)*cos(2*theta))*(sqrt(g0sq)*sqrt(g0sq)))/(2.*(g0sq*g0sq)*(Mh2*Mh2)) + 
        (3*a2*ct*A0(Mh1)*cos(3*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(8.*g0sq) + 
        (3*a2*ct*(Mh2*Mh2)*A0(Mh1)*cos(3*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(8.*g0sq*(Mh1*Mh1)) - 
        (3*ct*(Mh1*Mh1)*A0(Mh1)*cos(3*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(64.*(MW*MW)) - 
        (3*ct*(Mh2*Mh2)*A0(Mh1)*cos(3*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(32.*(MW*MW)) - 
        (3*(a2*a2)*ct*(MW*MW)*A0(Mh1)*cos(3*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(4.*(g0sq*g0sq)*(Mh1*Mh1)) - 
        (b3*st*A0(Mh1)*cos(3*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(16.*sqrt(g0sq)*MW) - 
        (b3*(Mh1*Mh1)*st*A0(Mh1)*cos(3*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(8.*sqrt(g0sq)*(Mh2*Mh2)*MW) + 
        (a2*b3*MW*st*A0(Mh1)*cos(3*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(4.*g0sq*sqrt(g0sq)*(Mh2*Mh2)) - 
        (3*(b3*b3)*ct*A0(Mh2)*cos(3*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(4.*g0sq*(Mh2*Mh2)) - 
        (3*b3*st*A0(Mh2)*cos(3*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(8.*sqrt(g0sq)*MW) - 
        (3*a2*b3*MW*st*A0(Mh2)*cos(3*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(2.*g0sq*sqrt(g0sq)*(Mh2*Mh2)) + 
        (3*b3*(Mt*Mt)*st*A0(Mt)*cos(3*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(sqrt(g0sq)*(Mh2*Mh2)*MW) - 
        (b3*st*A0(MW)*cos(3*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(4.*sqrt(g0sq)*MW) - 
        (b3*st*A0(MZ)*cos(3*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(8.*sqrt(g0sq)*MW) - 
        (3*a2*ct*A0(Mh1)*cos(2*theta)*cos(3*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(4.*g0sq) - 
        (3*a2*ct*(Mh2*Mh2)*A0(Mh1)*cos(2*theta)*cos(3*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(8.*g0sq*(Mh1*Mh1)) + 
        (3*ct*(Mh1*Mh1)*A0(Mh1)*cos(2*theta)*cos(3*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(64.*(MW*MW)) + 
        (3*ct*(Mh2*Mh2)*A0(Mh1)*cos(2*theta)*cos(3*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(32.*(MW*MW)) + 
        (9*(a2*a2)*ct*(MW*MW)*A0(Mh1)*cos(2*theta)*cos(3*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(4.*(g0sq*g0sq)*(Mh1*Mh1)) - 
        (b3*st*A0(Mh1)*cos(2*theta)*cos(3*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(16.*sqrt(g0sq)*MW) - 
        (b3*(Mh1*Mh1)*st*A0(Mh1)*cos(2*theta)*cos(3*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(8.*sqrt(g0sq)*(Mh2*Mh2)*MW) + 
        (3*a2*b3*MW*st*A0(Mh1)*cos(2*theta)*cos(3*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(4.*g0sq*sqrt(g0sq)*(Mh2*Mh2)) + 
        (3*b3*st*A0(Mh2)*cos(2*theta)*cos(3*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(8.*sqrt(g0sq)*MW) - 
        (3*a2*b3*MW*st*A0(Mh2)*cos(2*theta)*cos(3*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(2.*g0sq*sqrt(g0sq)*(Mh2*Mh2)) - 
        (b3*b3*A0(Mh2)*(cos(3*theta)*cos(3*theta))*(sqrt(g0sq)*sqrt(g0sq)))/(8.*g0sq*(Mh2*Mh2)) + 
        (9*b3*(ct*ct)*A0(Mh1)*sin(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(16.*sqrt(g0sq)*MW) + 
        (3*a2*b3*(ct*ct)*MW*A0(Mh1)*sin(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(4.*g0sq*sqrt(g0sq)*(Mh1*Mh1)) - 
        (3*(b3*b3)*ct*st*A0(Mh1)*sin(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(4.*g0sq*(Mh2*Mh2)) - 
        (3*b3*(st*st)*A0(Mh1)*sin(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(8.*sqrt(g0sq)*MW) - 
        (3*a2*b3*MW*(st*st)*A0(Mh1)*sin(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(2.*g0sq*sqrt(g0sq)*(Mh2*Mh2)) - 
        (b3*b3*ct*(st*st*st)*A0(Mh1)*sin(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(g0sq*(Mh1*Mh1)) + 
        (b3*(ct*ct)*A0(Mh2)*sin(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(4.*sqrt(g0sq)*MW) + 
        (b3*(ct*ct)*(Mh2*Mh2)*A0(Mh2)*sin(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(2.*sqrt(g0sq)*(Mh1*Mh1)*MW) - 
        (a2*b3*(ct*ct)*MW*A0(Mh2)*sin(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(g0sq*sqrt(g0sq)*(Mh1*Mh1)) - 
        (6*b3*(ct*ct)*(Mt*Mt)*A0(Mt)*sin(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(sqrt(g0sq)*(Mh1*Mh1)*MW) + 
        (b3*(ct*ct)*A0(MW)*sin(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(2.*sqrt(g0sq)*MW) + 
        (b3*(ct*ct)*A0(MZ)*sin(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(4.*sqrt(g0sq)*MW) + 
        (3*b3*(st*st)*A0(Mh1)*cos(2*theta)*sin(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(8.*sqrt(g0sq)*MW) - 
        (3*a2*b3*MW*(st*st)*A0(Mh1)*cos(2*theta)*sin(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(2.*g0sq*sqrt(g0sq)*(Mh2*Mh2)) - 
        (b3*(ct*ct)*A0(Mh2)*cos(2*theta)*sin(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(4.*sqrt(g0sq)*MW) - 
        (b3*(ct*ct)*(Mh2*Mh2)*A0(Mh2)*cos(2*theta)*sin(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(2.*sqrt(g0sq)*(Mh1*Mh1)*MW) + 
        (3*a2*b3*(ct*ct)*MW*A0(Mh2)*cos(2*theta)*sin(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(g0sq*sqrt(g0sq)*(Mh1*Mh1)) + 
        (3*b3*ct*A0(Mh1)*cos(3*theta)*sin(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(16.*sqrt(g0sq)*MW) - 
        (3*a2*b3*ct*MW*A0(Mh1)*cos(3*theta)*sin(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(4.*g0sq*sqrt(g0sq)*(Mh1*Mh1)) - 
        (b3*b3*st*A0(Mh1)*cos(3*theta)*sin(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(4.*g0sq*(Mh2*Mh2)) - 
        (b3*b3*(ct*ct)*A0(Mh2)*(sin(2*theta)*sin(2*theta))*(sqrt(g0sq)*sqrt(g0sq)))/(2.*g0sq*(Mh1*Mh1))
	);

	Float res = others.fin() + epsPart.div();
   	return (double) res;
}


double Renormalization::CalcSelfEnergyW(double k) {

	// Part proportional to epsilon
	Integral epsPart = ( 
        -2*g0sq*A0(MW) + (2*(ct*ct)*g0sq*(MW*MW)*A0(MW))/(Mh1*Mh1) + (2*g0sq*(MW*MW)*(st*st)*A0(MW))/(Mh2*Mh2) + 
		(ct*ct*g0sq*(MZ*MZ)*A0(MZ))/(Mh1*Mh1) - (2*(g0sq*g0sq)*A0(MZ))/(g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW)) + 
		(g0sq*(MZ*MZ)*(st*st)*A0(MZ))/(Mh2*Mh2) - (8*(g0sq*g0sq)*B00(k,MW,0))/(g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW)) + 
		(8*(g0sq*g0sq)*(MZ*MZ)*B00(k,MW,0))/(MW*MW*(g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW))) + 
		(8*(g0sq*g0sq)*B00(k,MW,MZ))/(g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW))
    );
	
	Integral others = (
		(-5*(ct*ct)*g0sq*A0(Mh1))/16. - (3*a2*(ct*ct)*(MW*MW)*A0(Mh1))/(4.*(Mh1*Mh1)) - (g0sq*(st*st)*A0(Mh1))/8. - 
		(g0sq*(Mh1*Mh1)*(st*st)*A0(Mh1))/(4.*(Mh2*Mh2)) + (a2*(MW*MW)*(st*st)*A0(Mh1))/(2.*(Mh2*Mh2)) + 
		(b3*ct*sqrt(g0sq)*MW*(st*st*st)*A0(Mh1))/(Mh1*Mh1) - (ct*ct*g0sq*A0(Mh2))/8. - 
		(ct*ct*g0sq*(Mh2*Mh2)*A0(Mh2))/(4.*(Mh1*Mh1)) + (a2*(ct*ct)*(MW*MW)*A0(Mh2))/(2.*(Mh1*Mh1)) - 
		(3*b3*ct*sqrt(g0sq)*MW*st*A0(Mh2))/(4.*(Mh2*Mh2)) - (g0sq*(st*st)*A0(Mh2))/8. - 
		(3*a2*(MW*MW)*(st*st)*A0(Mh2))/(2.*(Mh2*Mh2)) + (6*(ct*ct)*g0sq*(Mt*Mt)*A0(Mt))/(Mh1*Mh1) + 
		(6*g0sq*(Mt*Mt)*(st*st)*A0(Mt))/(Mh2*Mh2) + (7*g0sq*A0(MW))/2. - (ct*ct*g0sq*A0(MW))/2. - 
		(3*(ct*ct)*g0sq*(MW*MW)*A0(MW))/(Mh1*Mh1) - (g0sq*(st*st)*A0(MW))/2. - (3*g0sq*(MW*MW)*(st*st)*A0(MW))/(Mh2*Mh2) + 
		(g0sq*A0(MZ))/4. - (ct*ct*g0sq*A0(MZ))/4. - (3*(ct*ct)*g0sq*(MZ*MZ)*A0(MZ))/(2.*(Mh1*Mh1)) + 
		(3*(g0sq*g0sq)*A0(MZ))/(g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW)) - (g0sq*(st*st)*A0(MZ))/4. - 
		(3*g0sq*(MZ*MZ)*(st*st)*A0(MZ))/(2.*(Mh2*Mh2)) + ct*ct*g0sq*(MW*MW)*B0(k,Mh1,MW) + 
		g0sq*(MW*MW)*(st*st)*B0(k,Mh2,MW) - (g0sq*g0sq*(MW*MW)*B0(k,MW,0))/(g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW)) + 
		(g0sq*g0sq*(MZ*MZ)*B0(k,MW,0))/(g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW)) + 
		(g0sq*g0sq*(MW*MW)*B0(k,MW,MZ))/(g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW)) - 
		(2*(g0sq*g0sq)*(MZ*MZ)*B0(k,MW,MZ))/(g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW)) + 
		(g0sq*g0sq*(MZ*MZ*MZ*MZ)*B0(k,MW,MZ))/(MW*MW*(g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW))) + 18*g0sq*B00(k,0,0) + 
		6*g0sq*B00(k,0,Mt) - ct*ct*g0sq*B00(k,Mh1,MW) - g0sq*(st*st)*B00(k,Mh2,MW) + 
		(8*(g0sq*g0sq)*B00(k,MW,0))/(g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW)) - 
		(8*(g0sq*g0sq)*(MZ*MZ)*B00(k,MW,0))/(MW*MW*(g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW))) - 
		(8*(g0sq*g0sq)*B00(k,MW,MZ))/(g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW)) - g0sq*B00(k,MZ,MW) - 
		(g0sq*(st*st)*A0(Mh1)*cos(2*theta))/8. - (g0sq*(Mh1*Mh1)*(st*st)*A0(Mh1)*cos(2*theta))/(4.*(Mh2*Mh2)) + 
		(3*a2*(MW*MW)*(st*st)*A0(Mh1)*cos(2*theta))/(2.*(Mh2*Mh2)) + (ct*ct*g0sq*A0(Mh2)*cos(2*theta))/8. + 
		(ct*ct*g0sq*(Mh2*Mh2)*A0(Mh2)*cos(2*theta))/(4.*(Mh1*Mh1)) - 
		(3*a2*(ct*ct)*(MW*MW)*A0(Mh2)*cos(2*theta))/(2.*(Mh1*Mh1)) + (3*g0sq*(st*st)*A0(Mh2)*cos(2*theta))/8. - 
		(3*a2*(MW*MW)*(st*st)*A0(Mh2)*cos(2*theta))/(2.*(Mh2*Mh2)) - (3*ct*g0sq*A0(Mh1)*cos(3*theta))/16. + 
		(3*a2*ct*(MW*MW)*A0(Mh1)*cos(3*theta))/(4.*(Mh1*Mh1)) - (b3*sqrt(g0sq)*MW*st*A0(Mh2)*cos(3*theta))/(4.*(Mh2*Mh2)) - 
		9*g0sq*K(k,0,0,1,1,0) - 3*g0sq*K(k,0,Mt,1,1,0) + 
		(g0sq*g0sq*K(k,MW,0,2,2,5*(k*k)))/(g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW)) - 
		(g0sq*g0sq*(MZ*MZ)*K(k,MW,0,2,2,5*(k*k)))/(MW*MW*(g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW))) - 
		(g0sq*g0sq*K(k,MW,MZ,2,2,5*(k*k)))/(g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW)) - 
		(b3*sqrt(g0sq)*MW*(st*st)*A0(Mh1)*sin(2*theta))/(2.*(Mh2*Mh2)) + 
		(b3*(ct*ct)*sqrt(g0sq)*MW*A0(Mh2)*sin(2*theta))/(2.*(Mh1*Mh1))
   );

   Float res = others.fin() + epsPart.div();
   return (double) res;
}


double Renormalization::CalcSelfEnergyZ(double k) {

	// Part proportional to epsilon
	Integral epsPart = (2*(ct*ct)*g0sq*(MZ*MZ)*A0(MW))/(Mh1*Mh1) - (4*(g0sq*g0sq)*A0(MW))/(g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW)) + 
		(2*g0sq*(MZ*MZ)*(st*st)*A0(MW))/(Mh2*Mh2) + (ct*ct*g0sq*(MZ*MZ*MZ*MZ)*A0(MZ))/(Mh1*Mh1*(MW*MW)) + 
		(g0sq*(MZ*MZ*MZ*MZ)*(st*st)*A0(MZ))/(Mh2*Mh2*(MW*MW)) + 
		(8*(g0sq*g0sq)*B00(k,MW,MW))/(g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW));

	Integral others = (
		(ct*ct*g0sq*(MZ*MZ)*A0(Mh1))/(4.*(MW*MW)) + (g0sq*(MZ*MZ)*(st*st)*A0(Mh2))/(4.*(MW*MW)) - 
		(3*(ct*ct)*g0sq*(MZ*MZ)*A0(MW))/(Mh1*Mh1) + (8*(g0sq*g0sq)*A0(MW))/(g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW)) - 
		(2*(g0sq*g0sq)*(MZ*MZ)*A0(MW))/(MW*MW*(g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW))) + 
		(g0sq*g0sq*(MZ*MZ*MZ*MZ)*A0(MW))/(2.*(MW*MW*MW*MW)*(g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW))) - 
		(3*g0sq*(MZ*MZ)*(st*st)*A0(MW))/(Mh2*Mh2) + (g0sq*(MZ*MZ)*A0(MZ))/(4.*(MW*MW)) - 
		(3*(ct*ct)*g0sq*(MZ*MZ*MZ*MZ)*A0(MZ))/(2.*(Mh1*Mh1)*(MW*MW)) - 
		(3*g0sq*(MZ*MZ*MZ*MZ)*(st*st)*A0(MZ))/(2.*(Mh2*Mh2)*(MW*MW)) + (ct*ct*g0sq*(MZ*MZ*MZ*MZ)*B0(k,Mh1,MZ))/(MW*MW) + 
		(g0sq*(MZ*MZ*MZ*MZ)*(st*st)*B0(k,Mh2,MZ))/(MW*MW) + 
		(2*(g0sq*g0sq)*(MW*MW)*B0(k,MW,MW))/(g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW)) - 
		(4*(g0sq*g0sq)*(MZ*MZ)*B0(k,MW,MW))/(g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW)) + 
		(2*(g0sq*g0sq)*(MZ*MZ*MZ*MZ)*B0(k,MW,MW))/(MW*MW*(g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW))) + 
		(21*g0sq*(MZ*MZ)*B00(k,0,0))/(MW*MW) + (160*(g0sq*g0sq*g0sq)*(MZ*MZ)*B00(k,0,0))/
			(3.*(MW*MW)*((g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW))*(g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW)))) - 
		(320*(g0sq*g0sq*g0sq)*(MZ*MZ*MZ*MZ)*B00(k,0,0))/
			(3.*(MW*MW*MW*MW)*((g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW))*(g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW)))) + 
		(160*(g0sq*g0sq*g0sq)*(MZ*MZ*MZ*MZ*MZ*MZ)*B00(k,0,0))/
			(3.*(MW*MW*MW*MW*MW*MW)*((g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW))*(g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW)))) + 
		(40*(g0sq*g0sq)*(MZ*MZ)*B00(k,0,0))/(MW*MW*(g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW))) - 
		(40*(g0sq*g0sq)*(MZ*MZ*MZ*MZ)*B00(k,0,0))/(MW*MW*MW*MW*(g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW))) - 
		(ct*ct*g0sq*(MZ*MZ)*B00(k,Mh1,MZ))/(MW*MW) - (g0sq*(MZ*MZ)*(st*st)*B00(k,Mh2,MZ))/(MW*MW) + 
		(3*g0sq*(MZ*MZ)*B00(k,Mt,Mt))/(MW*MW) + (32*(g0sq*g0sq*g0sq)*(MZ*MZ)*B00(k,Mt,Mt))/
			(3.*(MW*MW)*((g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW))*(g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW)))) - 
		(64*(g0sq*g0sq*g0sq)*(MZ*MZ*MZ*MZ)*B00(k,Mt,Mt))/
			(3.*(MW*MW*MW*MW)*((g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW))*(g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW)))) + 
		(32*(g0sq*g0sq*g0sq)*(MZ*MZ*MZ*MZ*MZ*MZ)*B00(k,Mt,Mt))/
			(3.*(MW*MW*MW*MW*MW*MW)*((g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW))*(g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW)))) + 
		(8*(g0sq*g0sq)*(MZ*MZ)*B00(k,Mt,Mt))/(MW*MW*(g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW))) - 
		(8*(g0sq*g0sq)*(MZ*MZ*MZ*MZ)*B00(k,Mt,Mt))/(MW*MW*MW*MW*(g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW))) - 
		(12*(g0sq*g0sq)*B00(k,MW,MW))/(g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW)) + 
		(4*(g0sq*g0sq)*(MZ*MZ)*B00(k,MW,MW))/(MW*MW*(g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW))) - 
		(g0sq*g0sq*(MZ*MZ*MZ*MZ)*B00(k,MW,MW))/(MW*MW*MW*MW*(g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW))) - 
		(21*g0sq*(MZ*MZ)*K(k,0,0,1,1,0))/(2.*(MW*MW)) - 
		(80*(g0sq*g0sq*g0sq)*(MZ*MZ)*K(k,0,0,1,1,0))/
			(3.*(MW*MW)*((g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW))*(g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW)))) + 
		(160*(g0sq*g0sq*g0sq)*(MZ*MZ*MZ*MZ)*K(k,0,0,1,1,0))/
			(3.*(MW*MW*MW*MW)*((g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW))*(g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW)))) - 
		(80*(g0sq*g0sq*g0sq)*(MZ*MZ*MZ*MZ*MZ*MZ)*K(k,0,0,1,1,0))/
			(3.*(MW*MW*MW*MW*MW*MW)*((g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW))*(g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW)))) - 
		(20*(g0sq*g0sq)*(MZ*MZ)*K(k,0,0,1,1,0))/(MW*MW*(g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW))) + 
		(20*(g0sq*g0sq)*(MZ*MZ*MZ*MZ)*K(k,0,0,1,1,0))/(MW*MW*MW*MW*(g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW))) - 
		(3*g0sq*(MZ*MZ)*K(k,Mt,Mt,1,1,(-8*g0sq*(Mt*Mt)*(-(MW*MW) + MZ*MZ)*
				(-1 + (4*g0sq*(-(MW*MW) + MZ*MZ))/(3.*(MW*MW)*(g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW)))))/
				(3.*(MW*MW)*(g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW))*
				(1 + (8*g0sq*(-(MW*MW) + MZ*MZ)*(-1 + (4*g0sq*(-(MW*MW) + MZ*MZ))/
						(3.*(MW*MW)*(g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW)))))/
					(3.*(MW*MW)*(g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW)))))))/(2.*(MW*MW)) - 
		(16*(g0sq*g0sq*g0sq)*(MZ*MZ)*K(k,Mt,Mt,1,1,(-8*g0sq*(Mt*Mt)*(-(MW*MW) + MZ*MZ)*
				(-1 + (4*g0sq*(-(MW*MW) + MZ*MZ))/(3.*(MW*MW)*(g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW)))))/
				(3.*(MW*MW)*(g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW))*
				(1 + (8*g0sq*(-(MW*MW) + MZ*MZ)*(-1 + (4*g0sq*(-(MW*MW) + MZ*MZ))/
						(3.*(MW*MW)*(g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW)))))/
					(3.*(MW*MW)*(g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW)))))))/
			(3.*(MW*MW)*((g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW))*(g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW)))) + 
		(32*(g0sq*g0sq*g0sq)*(MZ*MZ*MZ*MZ)*K(k,Mt,Mt,1,1,
			(-8*g0sq*(Mt*Mt)*(-(MW*MW) + MZ*MZ)*(-1 + 
					(4*g0sq*(-(MW*MW) + MZ*MZ))/(3.*(MW*MW)*(g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW)))))/
				(3.*(MW*MW)*(g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW))*
				(1 + (8*g0sq*(-(MW*MW) + MZ*MZ)*(-1 + (4*g0sq*(-(MW*MW) + MZ*MZ))/
						(3.*(MW*MW)*(g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW)))))/
					(3.*(MW*MW)*(g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW)))))))/
			(3.*(MW*MW*MW*MW)*((g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW))*(g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW)))) - 
		(16*(g0sq*g0sq*g0sq)*(MZ*MZ*MZ*MZ*MZ*MZ)*K(k,Mt,Mt,1,1,
			(-8*g0sq*(Mt*Mt)*(-(MW*MW) + MZ*MZ)*(-1 + 
					(4*g0sq*(-(MW*MW) + MZ*MZ))/(3.*(MW*MW)*(g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW)))))/
				(3.*(MW*MW)*(g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW))*
				(1 + (8*g0sq*(-(MW*MW) + MZ*MZ)*(-1 + (4*g0sq*(-(MW*MW) + MZ*MZ))/
						(3.*(MW*MW)*(g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW)))))/
					(3.*(MW*MW)*(g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW)))))))/
			(3.*(MW*MW*MW*MW*MW*MW)*((g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW))*(g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW)))) - 
		(4*(g0sq*g0sq)*(MZ*MZ)*K(k,Mt,Mt,1,1,(-8*g0sq*(Mt*Mt)*(-(MW*MW) + MZ*MZ)*
				(-1 + (4*g0sq*(-(MW*MW) + MZ*MZ))/(3.*(MW*MW)*(g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW)))))/
				(3.*(MW*MW)*(g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW))*
				(1 + (8*g0sq*(-(MW*MW) + MZ*MZ)*(-1 + (4*g0sq*(-(MW*MW) + MZ*MZ))/
						(3.*(MW*MW)*(g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW)))))/
					(3.*(MW*MW)*(g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW)))))))/
			(MW*MW*(g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW))) + 
		(4*(g0sq*g0sq)*(MZ*MZ*MZ*MZ)*K(k,Mt,Mt,1,1,(-8*g0sq*(Mt*Mt)*(-(MW*MW) + MZ*MZ)*
				(-1 + (4*g0sq*(-(MW*MW) + MZ*MZ))/(3.*(MW*MW)*(g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW)))))/
				(3.*(MW*MW)*(g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW))*
				(1 + (8*g0sq*(-(MW*MW) + MZ*MZ)*(-1 + (4*g0sq*(-(MW*MW) + MZ*MZ))/
						(3.*(MW*MW)*(g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW)))))/
					(3.*(MW*MW)*(g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW)))))))/
			(MW*MW*MW*MW*(g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW))) - 
		(g0sq*g0sq*K(k,MW,MW,2,2,5*(k*k)))/(g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW)) - 
		(3*a2*(ct*ct)*(MZ*MZ)*A0(Mh1)*sqrt(g0sq))/(4.*sqrt(g0sq)*(Mh1*Mh1)) - 
		(9*(ct*ct)*sqrt(g0sq)*(MZ*MZ)*A0(Mh1)*sqrt(g0sq))/(16.*(MW*MW)) + 
		(a2*(MZ*MZ)*(st*st)*A0(Mh1)*sqrt(g0sq))/(2.*sqrt(g0sq)*(Mh2*Mh2)) - 
		(sqrt(g0sq)*(MZ*MZ)*(st*st)*A0(Mh1)*sqrt(g0sq))/(8.*(MW*MW)) - 
		(sqrt(g0sq)*(Mh1*Mh1)*(MZ*MZ)*(st*st)*A0(Mh1)*sqrt(g0sq))/(4.*(Mh2*Mh2)*(MW*MW)) + 
		(b3*ct*(MZ*MZ)*(st*st*st)*A0(Mh1)*sqrt(g0sq))/(Mh1*Mh1*MW) + 
		(a2*(ct*ct)*(MZ*MZ)*A0(Mh2)*sqrt(g0sq))/(2.*sqrt(g0sq)*(Mh1*Mh1)) - 
		(ct*ct*sqrt(g0sq)*(MZ*MZ)*A0(Mh2)*sqrt(g0sq))/(8.*(MW*MW)) - 
		(ct*ct*sqrt(g0sq)*(Mh2*Mh2)*(MZ*MZ)*A0(Mh2)*sqrt(g0sq))/(4.*(Mh1*Mh1)*(MW*MW)) - 
		(3*b3*ct*(MZ*MZ)*st*A0(Mh2)*sqrt(g0sq))/(4.*(Mh2*Mh2)*MW) - 
		(3*a2*(MZ*MZ)*(st*st)*A0(Mh2)*sqrt(g0sq))/(2.*sqrt(g0sq)*(Mh2*Mh2)) - 
		(3*sqrt(g0sq)*(MZ*MZ)*(st*st)*A0(Mh2)*sqrt(g0sq))/(8.*(MW*MW)) + 
		(6*(ct*ct)*sqrt(g0sq)*(Mt*Mt)*(MZ*MZ)*A0(Mt)*sqrt(g0sq))/(Mh1*Mh1*(MW*MW)) + 
		(6*sqrt(g0sq)*(Mt*Mt)*(MZ*MZ)*(st*st)*A0(Mt)*sqrt(g0sq))/(Mh2*Mh2*(MW*MW)) - 
		(ct*ct*sqrt(g0sq)*(MZ*MZ)*A0(MW)*sqrt(g0sq))/(2.*(MW*MW)) - 
		(sqrt(g0sq)*(MZ*MZ)*(st*st)*A0(MW)*sqrt(g0sq))/(2.*(MW*MW)) - 
		(ct*ct*sqrt(g0sq)*(MZ*MZ)*A0(MZ)*sqrt(g0sq))/(4.*(MW*MW)) - 
		(sqrt(g0sq)*(MZ*MZ)*(st*st)*A0(MZ)*sqrt(g0sq))/(4.*(MW*MW)) + 
		(3*a2*(MZ*MZ)*(st*st)*A0(Mh1)*cos(2*theta)*sqrt(g0sq))/(2.*sqrt(g0sq)*(Mh2*Mh2)) - 
		(sqrt(g0sq)*(MZ*MZ)*(st*st)*A0(Mh1)*cos(2*theta)*sqrt(g0sq))/(8.*(MW*MW)) - 
		(sqrt(g0sq)*(Mh1*Mh1)*(MZ*MZ)*(st*st)*A0(Mh1)*cos(2*theta)*sqrt(g0sq))/(4.*(Mh2*Mh2)*(MW*MW)) - 
		(3*a2*(ct*ct)*(MZ*MZ)*A0(Mh2)*cos(2*theta)*sqrt(g0sq))/(2.*sqrt(g0sq)*(Mh1*Mh1)) + 
		(ct*ct*sqrt(g0sq)*(MZ*MZ)*A0(Mh2)*cos(2*theta)*sqrt(g0sq))/(8.*(MW*MW)) + 
		(ct*ct*sqrt(g0sq)*(Mh2*Mh2)*(MZ*MZ)*A0(Mh2)*cos(2*theta)*sqrt(g0sq))/(4.*(Mh1*Mh1)*(MW*MW)) - 
		(3*a2*(MZ*MZ)*(st*st)*A0(Mh2)*cos(2*theta)*sqrt(g0sq))/(2.*sqrt(g0sq)*(Mh2*Mh2)) + 
		(3*sqrt(g0sq)*(MZ*MZ)*(st*st)*A0(Mh2)*cos(2*theta)*sqrt(g0sq))/(8.*(MW*MW)) + 
		(3*a2*ct*(MZ*MZ)*A0(Mh1)*cos(3*theta)*sqrt(g0sq))/(4.*sqrt(g0sq)*(Mh1*Mh1)) - 
		(3*ct*sqrt(g0sq)*(MZ*MZ)*A0(Mh1)*cos(3*theta)*sqrt(g0sq))/(16.*(MW*MW)) - 
		(b3*(MZ*MZ)*st*A0(Mh2)*cos(3*theta)*sqrt(g0sq))/(4.*(Mh2*Mh2)*MW) - 
		(b3*(MZ*MZ)*(st*st)*A0(Mh1)*sin(2*theta)*sqrt(g0sq))/(2.*(Mh2*Mh2)*MW) + 
		(b3*(ct*ct)*(MZ*MZ)*A0(Mh2)*sin(2*theta)*sqrt(g0sq))/(2.*(Mh1*Mh1)*MW)
	);

    Float res = others.fin() + epsPart.div();
  	return (double) res;
}

double Renormalization::CalcSelfEnergyTop(double k) {

    // Calculate scalar + vector part of top quark self energy

	// Part proportional to epsilon
	Integral epsPart = (
        (8*gs2*B0(k,Mt,0))/3. - (8*(g0sq*g0sq)*B0(k,Mt,0))/(9.*(g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW))) + 
        (8*(g0sq*g0sq)*(MZ*MZ)*B0(k,Mt,0))/(9.*(MW*MW)*(g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW))) + 
        (8*(g0sq*g0sq)*B0(k,Mt,MZ))/(9.*(g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW))) - 
        (10*(g0sq*g0sq)*(MZ*MZ)*B0(k,Mt,MZ))/(9.*(MW*MW)*(g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW))) + 
        (2*(g0sq*g0sq)*(MZ*MZ*MZ*MZ)*B0(k,Mt,MZ))/(9.*(MW*MW*MW*MW)*(g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW))) + 
        (g0sq*B1(k,0,MW))/2. + (8*gs2*B1(k,Mt,0))/3. - 
        (8*(g0sq*g0sq)*B1(k,Mt,0))/(9.*(g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW))) + 
        (8*(g0sq*g0sq)*(MZ*MZ)*B1(k,Mt,0))/(9.*(MW*MW)*(g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW))) + 
        (8*(g0sq*g0sq)*B1(k,Mt,MZ))/(9.*(g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW))) - 
        (10*(g0sq*g0sq)*(MZ*MZ)*B1(k,Mt,MZ))/(9.*(MW*MW)*(g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW))) + 
        (17*(g0sq*g0sq)*(MZ*MZ*MZ*MZ)*B1(k,Mt,MZ))/(36.*(MW*MW*MW*MW)*(g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW))) + 
        (ct*ct*sqrt(g0sq)*A0(MW)*sqrt(g0sq))/(Mh1*Mh1) + (sqrt(g0sq)*(st*st)*A0(MW)*sqrt(g0sq))/(Mh2*Mh2) + 
        (ct*ct*sqrt(g0sq)*(MZ*MZ)*A0(MZ)*sqrt(g0sq))/(2.*(Mh1*Mh1)*(MW*MW)) + 
        (sqrt(g0sq)*(MZ*MZ)*(st*st)*A0(MZ)*sqrt(g0sq))/(2.*(Mh2*Mh2)*(MW*MW))
    );

    Integral others = ( 
        (-16*gs2*B0(k,Mt,0))/3. + (16*(g0sq*g0sq)*B0(k,Mt,0))/(9.*(g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW))) - 
        (16*(g0sq*g0sq)*(MZ*MZ)*B0(k,Mt,0))/(9.*(MW*MW)*(g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW))) + 
        (ct*ct*g0sq*(Mt*Mt)*B0(k,Mt,Mh1))/(4.*(MW*MW)) + (g0sq*(Mt*Mt)*(st*st)*B0(k,Mt,Mh2))/(4.*(MW*MW)) - 
        (g0sq*(Mt*Mt)*B0(k,Mt,MZ))/(4.*(MW*MW)) - (16*(g0sq*g0sq)*B0(k,Mt,MZ))/
        (9.*(g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW))) + 
        (20*(g0sq*g0sq)*(MZ*MZ)*B0(k,Mt,MZ))/(9.*(MW*MW)*(g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW))) - 
        (4*(g0sq*g0sq)*(MZ*MZ*MZ*MZ)*B0(k,Mt,MZ))/(9.*(MW*MW*MW*MW)*(g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW))) - 
        (g0sq*B1(k,0,MW))/2. - (g0sq*(Mt*Mt)*B1(k,0,MW))/(4.*(MW*MW)) - (8*gs2*B1(k,Mt,0))/3. + 
        (8*(g0sq*g0sq)*B1(k,Mt,0))/(9.*(g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW))) - 
        (8*(g0sq*g0sq)*(MZ*MZ)*B1(k,Mt,0))/(9.*(MW*MW)*(g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW))) - 
        (ct*ct*g0sq*(Mt*Mt)*B1(k,Mt,Mh1))/(4.*(MW*MW)) - (g0sq*(Mt*Mt)*(st*st)*B1(k,Mt,Mh2))/(4.*(MW*MW)) - 
        (g0sq*(Mt*Mt)*B1(k,Mt,MZ))/(4.*(MW*MW)) - (8*(g0sq*g0sq)*B1(k,Mt,MZ))/
        (9.*(g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW))) + 
        (10*(g0sq*g0sq)*(MZ*MZ)*B1(k,Mt,MZ))/(9.*(MW*MW)*(g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW))) - 
        (17*(g0sq*g0sq)*(MZ*MZ*MZ*MZ)*B1(k,Mt,MZ))/(36.*(MW*MW*MW*MW)*(g0sq + (g0sq*(-(MW*MW) + MZ*MZ))/(MW*MW))) - 
        (3*(ct*ct)*sqrt(g0sq)*A0(MW)*sqrt(g0sq))/(2.*(Mh1*Mh1)) - (3*sqrt(g0sq)*(st*st)*A0(MW)*sqrt(g0sq))/(2.*(Mh2*Mh2)) - 
        (3*(ct*ct)*sqrt(g0sq)*(MZ*MZ)*A0(MZ)*sqrt(g0sq))/(4.*(Mh1*Mh1)*(MW*MW)) - 
        (3*sqrt(g0sq)*(MZ*MZ)*(st*st)*A0(MZ)*sqrt(g0sq))/(4.*(Mh2*Mh2)*(MW*MW)) - 
        (3*a2*(ct*ct)*A0(Mh1)*(sqrt(g0sq)*sqrt(g0sq)))/(8.*g0sq*(Mh1*Mh1)) - 
        (9*(ct*ct)*A0(Mh1)*(sqrt(g0sq)*sqrt(g0sq)))/(32.*(MW*MW)) + 
        (a2*(st*st)*A0(Mh1)*(sqrt(g0sq)*sqrt(g0sq)))/(4.*g0sq*(Mh2*Mh2)) - 
        (st*st*A0(Mh1)*(sqrt(g0sq)*sqrt(g0sq)))/(16.*(MW*MW)) - 
        (Mh1*Mh1*(st*st)*A0(Mh1)*(sqrt(g0sq)*sqrt(g0sq)))/(8.*(Mh2*Mh2)*(MW*MW)) + 
        (b3*ct*(st*st*st)*A0(Mh1)*(sqrt(g0sq)*sqrt(g0sq)))/(2.*sqrt(g0sq)*(Mh1*Mh1)*MW) + 
        (a2*(ct*ct)*A0(Mh2)*(sqrt(g0sq)*sqrt(g0sq)))/(4.*g0sq*(Mh1*Mh1)) - 
        (ct*ct*A0(Mh2)*(sqrt(g0sq)*sqrt(g0sq)))/(16.*(MW*MW)) - 
        (ct*ct*(Mh2*Mh2)*A0(Mh2)*(sqrt(g0sq)*sqrt(g0sq)))/(8.*(Mh1*Mh1)*(MW*MW)) - 
        (3*b3*ct*st*A0(Mh2)*(sqrt(g0sq)*sqrt(g0sq)))/(8.*sqrt(g0sq)*(Mh2*Mh2)*MW) - 
        (3*a2*(st*st)*A0(Mh2)*(sqrt(g0sq)*sqrt(g0sq)))/(4.*g0sq*(Mh2*Mh2)) - 
        (3*(st*st)*A0(Mh2)*(sqrt(g0sq)*sqrt(g0sq)))/(16.*(MW*MW)) + 
        (3*(ct*ct)*(Mt*Mt)*A0(Mt)*(sqrt(g0sq)*sqrt(g0sq)))/(Mh1*Mh1*(MW*MW)) + 
        (3*(Mt*Mt)*(st*st)*A0(Mt)*(sqrt(g0sq)*sqrt(g0sq)))/(Mh2*Mh2*(MW*MW)) - 
        (ct*ct*A0(MW)*(sqrt(g0sq)*sqrt(g0sq)))/(4.*(MW*MW)) - (st*st*A0(MW)*(sqrt(g0sq)*sqrt(g0sq)))/(4.*(MW*MW)) - 
        (ct*ct*A0(MZ)*(sqrt(g0sq)*sqrt(g0sq)))/(8.*(MW*MW)) - (st*st*A0(MZ)*(sqrt(g0sq)*sqrt(g0sq)))/(8.*(MW*MW)) + 
        (3*a2*(st*st)*A0(Mh1)*cos(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(4.*g0sq*(Mh2*Mh2)) - 
        (st*st*A0(Mh1)*cos(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(16.*(MW*MW)) - 
        (Mh1*Mh1*(st*st)*A0(Mh1)*cos(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(8.*(Mh2*Mh2)*(MW*MW)) - 
        (3*a2*(ct*ct)*A0(Mh2)*cos(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(4.*g0sq*(Mh1*Mh1)) + 
        (ct*ct*A0(Mh2)*cos(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(16.*(MW*MW)) + 
        (ct*ct*(Mh2*Mh2)*A0(Mh2)*cos(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(8.*(Mh1*Mh1)*(MW*MW)) - 
        (3*a2*(st*st)*A0(Mh2)*cos(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(4.*g0sq*(Mh2*Mh2)) + 
        (3*(st*st)*A0(Mh2)*cos(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(16.*(MW*MW)) + 
        (3*a2*ct*A0(Mh1)*cos(3*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(8.*g0sq*(Mh1*Mh1)) - 
        (3*ct*A0(Mh1)*cos(3*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(32.*(MW*MW)) - 
        (b3*st*A0(Mh2)*cos(3*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(8.*sqrt(g0sq)*(Mh2*Mh2)*MW) - 
        (b3*(st*st)*A0(Mh1)*sin(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(4.*sqrt(g0sq)*(Mh2*Mh2)*MW) + 
        (b3*(ct*ct)*A0(Mh2)*sin(2*theta)*(sqrt(g0sq)*sqrt(g0sq)))/(4.*sqrt(g0sq)*(Mh1*Mh1)*MW)
    );

    Float res = others.fin() + epsPart.div();
  	return (double) res;
}