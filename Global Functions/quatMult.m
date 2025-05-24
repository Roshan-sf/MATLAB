function out = quatMult(p, q)
np = p(4);
nq = q(4);
epsP = p(1:3);
epsQ = q(1:3);

out = [np*epsQ + nq*epsP + vcross(epsP)*epsQ; np*nq - epsP'*epsQ];

end