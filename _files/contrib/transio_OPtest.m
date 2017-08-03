% transio_OPtest  - test transio overloaded operators
function transio_OPtest

% get three transio objects, two the same size, one 1x1 size
v1 = xff([neuroelf_path '/_files/colin/colin_brain_ICBMnorm.vmr'], 't');
t1 = v1.VMRData;
v1.ClearObject;
t1 = transio(filename(t1), 'le', 'uint8', offset(t1) + 256 * 256 * 124, [256, 256, 8]);
v2 = xff([neuroelf_path '/_files/colin/colin_brain.vmr'], 't');
t2 = v2.VMRData;
v2.ClearObject;
t2 = transio(filename(t2), 'le', 'uint8', offset(t2) + 256 * 256 * 124, [256, 256, 8]);
t3 = transio(filename(t1), 'le', 'uint8', offset(t1), [1, 1]);
clear v1 v2;

% perform operations on transio data
% EQ(OBJ,OBJ)
eq11 = t1 == t1;
eq12 = t1 == t2;
eq13 = t1 == t3;
eq21 = t2 == t1;
eq22 = t2 == t2;
eq23 = t2 == t3;
eq31 = t3 == t1;
eq32 = t3 == t2;
eq33 = t3 == t3;

% NE(OBJ,OBJ)
ne11 = t1 ~= t1;
ne12 = t1 ~= t2;
ne13 = t1 ~= t3;
ne21 = t2 ~= t1;
ne22 = t2 ~= t2;
ne23 = t2 ~= t3;
ne31 = t3 ~= t1;
ne32 = t3 ~= t2;
ne33 = t3 ~= t3;

% LT(OBJ,OBJ)
lt11 = t1 < t1;
lt12 = t1 < t2;
lt13 = t1 < t3;
lt21 = t2 < t1;
lt22 = t2 < t2;
lt23 = t2 < t3;
lt31 = t3 < t1;
lt32 = t3 < t2;
lt33 = t3 < t3;

% LE(OBJ,OBJ)
le11 = t1 <= t1;
le12 = t1 <= t2;
le13 = t1 <= t3;
le21 = t2 <= t1;
le22 = t2 <= t2;
le23 = t2 <= t3;
le31 = t3 <= t1;
le32 = t3 <= t2;
le33 = t3 <= t3;

% GE(OBJ,OBJ)
ge11 = t1 >= t1;
ge12 = t1 >= t2;
ge13 = t1 >= t3;
ge21 = t2 >= t1;
ge22 = t2 >= t2;
ge23 = t2 >= t3;
ge31 = t3 >= t1;
ge32 = t3 >= t2;
ge33 = t3 >= t3;

% GT(OBJ,OBJ)
gt11 = t1 > t1;
gt12 = t1 > t2;
gt13 = t1 > t3;
gt21 = t2 > t1;
gt22 = t2 > t2;
gt23 = t2 > t3;
gt31 = t3 > t1;
gt32 = t3 > t2;
gt33 = t3 > t3;

% perform operations with single value
eq1a = t1 == 0;
eq1b = 0 == t1;
eq3a = t3 == 0;
eq3b = 0 == t3;
ne1a = t1 ~= 0;
ne1b = 0 ~= t1;
ne3a = t3 ~= 0;
ne3b = 0 ~= t3;
lt1a = t1 < 0;
lt1b = 0 < t1;
lt3a = t3 < 0;
lt3b = 0 < t3;
le1a = t1 <= 0;
le1b = 0 <= t1;
le3a = t3 <= 0;
le3b = 0 <= t3;
ge1a = t1 >= 0;
ge1b = 0 >= t1;
ge3a = t3 >= 0;
ge3b = 0 >= t3;
gt1a = t1 > 0;
gt1b = 0 > t1;
gt3a = t3 > 0;
gt3b = 0 > t3;

% resolve to data
d1 = resolve(t1);
d2 = resolve(t2);
d3 = resolve(t3);

% EQ(OBJ,V) + EQ(V,OBJ)
eq11a = t1 == d1;
eq11b = d1 == t1;
eq12a = t1 == d2;
eq12b = d1 == t2;
eq13a = t1 == d3;
eq13b = d1 == t3;
eq21a = t2 == d1;
eq21b = d2 == t1;
eq22a = t2 == d2;
eq22b = d2 == t2;
eq23a = t2 == d3;
eq23b = d2 == t3;
eq31a = t3 == d1;
eq31b = d3 == t1;
eq32a = t3 == d2;
eq32b = d3 == t2;
eq33a = t3 == d3;
eq33b = d3 == t3;

% NE(OBJ,V) + NE(V,OBJ)
ne11a = t1 ~= d1;
ne11b = d1 ~= t1;
ne12a = t1 ~= d2;
ne12b = d1 ~= t2;
ne13a = t1 ~= d3;
ne13b = d1 ~= t3;
ne21a = t2 ~= d1;
ne21b = d2 ~= t1;
ne22a = t2 ~= d2;
ne22b = d2 ~= t2;
ne23a = t2 ~= d3;
ne23b = d2 ~= t3;
ne31a = t3 ~= d1;
ne31b = d3 ~= t1;
ne32a = t3 ~= d2;
ne32b = d3 ~= t2;
ne33a = t3 ~= d3;
ne33b = d3 ~= t3;

% LT(OBJ,V) + LT(V,OBJ)
lt11a = t1 < d1;
lt11b = d1 < t1;
lt12a = t1 < d2;
lt12b = d1 < t2;
lt13a = t1 < d3;
lt13b = d1 < t3;
lt21a = t2 < d1;
lt21b = d2 < t1;
lt22a = t2 < d2;
lt22b = d2 < t2;
lt23a = t2 < d3;
lt23b = d2 < t3;
lt31a = t3 < d1;
lt31b = d3 < t1;
lt32a = t3 < d2;
lt32b = d3 < t2;
lt33a = t3 < d3;
lt33b = d3 < t3;

% LE(OBJ,V) + LE(V,OBJ)
le11a = t1 <= d1;
le11b = d1 <= t1;
le12a = t1 <= d2;
le12b = d1 <= t2;
le13a = t1 <= d3;
le13b = d1 <= t3;
le21a = t2 <= d1;
le21b = d2 <= t1;
le22a = t2 <= d2;
le22b = d2 <= t2;
le23a = t2 <= d3;
le23b = d2 <= t3;
le31a = t3 <= d1;
le31b = d3 <= t1;
le32a = t3 <= d2;
le32b = d3 <= t2;
le33a = t3 <= d3;
le33b = d3 <= t3;

% LE(OBJ,V) + LE(V,OBJ)
ge11a = t1 >= d1;
ge11b = d1 >= t1;
ge12a = t1 >= d2;
ge12b = d1 >= t2;
ge13a = t1 >= d3;
ge13b = d1 >= t3;
ge21a = t2 >= d1;
ge21b = d2 >= t1;
ge22a = t2 >= d2;
ge22b = d2 >= t2;
ge23a = t2 >= d3;
ge23b = d2 >= t3;
ge31a = t3 >= d1;
ge31b = d3 >= t1;
ge32a = t3 >= d2;
ge32b = d3 >= t2;
ge33a = t3 >= d3;
ge33b = d3 >= t3;

% LT(OBJ,V) + LT(V,OBJ)
gt11a = t1 > d1;
gt11b = d1 > t1;
gt12a = t1 > d2;
gt12b = d1 > t2;
gt13a = t1 > d3;
gt13b = d1 > t3;
gt21a = t2 > d1;
gt21b = d2 > t1;
gt22a = t2 > d2;
gt22b = d2 > t2;
gt23a = t2 > d3;
gt23b = d2 > t3;
gt31a = t3 > d1;
gt31b = d3 > t1;
gt32a = t3 > d2;
gt32b = d3 > t2;
gt33a = t3 > d3;
gt33b = d3 > t3;

% perform the same operations on values
deq11 = d1 == d1;
deq12 = d1 == d2;
deq13 = d1 == d3;
deq21 = d2 == d1;
deq22 = d2 == d2;
deq23 = d2 == d3;
deq31 = d3 == d1;
deq32 = d3 == d2;
deq33 = d3 == d3;
dlt11 = d1 < d1;
dlt12 = d1 < d2;
dlt13 = d1 < d3;
dlt21 = d2 < d1;
dlt22 = d2 < d2;
dlt23 = d2 < d3;
dlt31 = d3 < d1;
dlt32 = d3 < d2;
dlt33 = d3 < d3;
dle11 = d1 <= d1;
dle12 = d1 <= d2;
dle13 = d1 <= d3;
dle21 = d2 <= d1;
dle22 = d2 <= d2;
dle23 = d2 <= d3;
dle31 = d3 <= d1;
dle32 = d3 <= d2;
dle33 = d3 <= d3;
dge11 = d1 >= d1;
dge12 = d1 >= d2;
dge13 = d1 >= d3;
dge21 = d2 >= d1;
dge22 = d2 >= d2;
dge23 = d2 >= d3;
dge31 = d3 >= d1;
dge32 = d3 >= d2;
dge33 = d3 >= d3;
dgt11 = d1 > d1;
dgt12 = d1 > d2;
dgt13 = d1 > d3;
dgt21 = d2 > d1;
dgt22 = d2 > d2;
dgt23 = d2 > d3;
dgt31 = d3 > d1;
dgt32 = d3 > d2;
dgt33 = d3 > d3;
dne11 = d1 ~= d1;
dne12 = d1 ~= d2;
dne13 = d1 ~= d3;
dne21 = d2 ~= d1;
dne22 = d2 ~= d2;
dne23 = d2 ~= d3;
dne31 = d3 ~= d1;
dne32 = d3 ~= d2;
dne33 = d3 ~= d3;
deq1a = d1 == 0;
deq1b = 0 == d1;
deq3a = d3 == 0;
deq3b = 0 == d3;
dne1a = d1 ~= 0;
dne1b = 0 ~= d1;
dne3a = d3 ~= 0;
dne3b = 0 ~= d3;
dge1a = d1 >= 0;
dge1b = 0 >= d1;
dge3a = d3 >= 0;
dge3b = 0 >= d3;
dgt1a = d1 > 0;
dgt1b = 0 > d1;
dgt3a = d3 > 0;
dgt3b = 0 > d3;
dlt1a = d1 < 0;
dlt1b = 0 < d1;
dlt3a = d3 < 0;
dlt3b = 0 < d3;
dle1a = d1 <= 0;
dle1b = 0 <= d1;
dle3a = d3 <= 0;
dle3b = 0 <= d3;
deq11a = d1 == d1;
deq11b = d1 == d1;
deq12a = d1 == d2;
deq12b = d1 == d2;
deq13a = d1 == d3;
deq13b = d1 == d3;
deq21a = d2 == d1;
deq21b = d2 == d1;
deq22a = d2 == d2;
deq22b = d2 == d2;
deq23a = d2 == d3;
deq23b = d2 == d3;
deq31a = d3 == d1;
deq31b = d3 == d1;
deq32a = d3 == d2;
deq32b = d3 == d2;
deq33a = d3 == d3;
deq33b = d3 == d3;
dne11a = d1 ~= d1;
dne11b = d1 ~= d1;
dne12a = d1 ~= d2;
dne12b = d1 ~= d2;
dne13a = d1 ~= d3;
dne13b = d1 ~= d3;
dne21a = d2 ~= d1;
dne21b = d2 ~= d1;
dne22a = d2 ~= d2;
dne22b = d2 ~= d2;
dne23a = d2 ~= d3;
dne23b = d2 ~= d3;
dne31a = d3 ~= d1;
dne31b = d3 ~= d1;
dne32a = d3 ~= d2;
dne32b = d3 ~= d2;
dne33a = d3 ~= d3;
dne33b = d3 ~= d3;
dlt11a = d1 < d1;
dlt11b = d1 < d1;
dlt12a = d1 < d2;
dlt12b = d1 < d2;
dlt13a = d1 < d3;
dlt13b = d1 < d3;
dlt21a = d2 < d1;
dlt21b = d2 < d1;
dlt22a = d2 < d2;
dlt22b = d2 < d2;
dlt23a = d2 < d3;
dlt23b = d2 < d3;
dlt31a = d3 < d1;
dlt31b = d3 < d1;
dlt32a = d3 < d2;
dlt32b = d3 < d2;
dlt33a = d3 < d3;
dlt33b = d3 < d3;
dle11a = d1 <= d1;
dle11b = d1 <= d1;
dle12a = d1 <= d2;
dle12b = d1 <= d2;
dle13a = d1 <= d3;
dle13b = d1 <= d3;
dle21a = d2 <= d1;
dle21b = d2 <= d1;
dle22a = d2 <= d2;
dle22b = d2 <= d2;
dle23a = d2 <= d3;
dle23b = d2 <= d3;
dle31a = d3 <= d1;
dle31b = d3 <= d1;
dle32a = d3 <= d2;
dle32b = d3 <= d2;
dle33a = d3 <= d3;
dle33b = d3 <= d3;
dge11a = d1 >= d1;
dge11b = d1 >= d1;
dge12a = d1 >= d2;
dge12b = d1 >= d2;
dge13a = d1 >= d3;
dge13b = d1 >= d3;
dge21a = d2 >= d1;
dge21b = d2 >= d1;
dge22a = d2 >= d2;
dge22b = d2 >= d2;
dge23a = d2 >= d3;
dge23b = d2 >= d3;
dge31a = d3 >= d1;
dge31b = d3 >= d1;
dge32a = d3 >= d2;
dge32b = d3 >= d2;
dge33a = d3 >= d3;
dge33b = d3 >= d3;
dgt11a = d1 > d1;
dgt11b = d1 > d1;
dgt12a = d1 > d2;
dgt12b = d1 > d2;
dgt13a = d1 > d3;
dgt13b = d1 > d3;
dgt21a = d2 > d1;
dgt21b = d2 > d1;
dgt22a = d2 > d2;
dgt22b = d2 > d2;
dgt23a = d2 > d3;
dgt23b = d2 > d3;
dgt31a = d3 > d1;
dgt31b = d3 > d1;
dgt32a = d3 > d2;
dgt32b = d3 > d2;
dgt33a = d3 > d3;
dgt33b = d3 > d3;

% clear actual values
clear d1 d2 d3;

% get list of names (starting with d)
w = whos;
wn = {w.name};
wn(cellfun('isempty', regexp(wn, '^d'))) = [];

% and compare results
for c = 1:numel(wn)
    eval(['if ~isequal(' wn{c} ', ' wn{c}(2:end) '), disp(wn{c}); end']);
end
