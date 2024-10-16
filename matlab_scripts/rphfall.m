%compute K(k) for Ray Pierrehumbert approach
k0_s = 0.5:0.5:10;
keff_s=[];
for k0=k0_s
  k0,fflush(1);
  rphf
  fflush(1);
  keff_s=[keff_s,kcoeff];
end

plot(0.5:0.5:10,keff_s)
xlabel("k");
titleb("Ke");
