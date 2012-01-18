function test_quadrature 


% Test LevelSymmetric

q = LevelSymmetric(8);
[mu eta] = angles(q);
[mu1 eta1] = angle(q, 26)
number_angles(q)


% Test UniformEqualWeight

q = UniformEqualWeight(8);
[mu eta] = angles(q);
[mu1 eta1] = angle(q, 26)
number_angles(q)


% Test QuadrupleRange

q = QuadrupleRange(8);
[mu eta] = angles(q);
[mu1 eta1] = angle(q, 26)
number_angles(q)







end