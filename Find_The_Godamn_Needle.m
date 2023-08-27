%% Find the Goddamn Needle
CompareSecondPartS = reshape(Ms*khatrirao(S, S)' + S*khatrirao(S, Ms)' + S*khatrirao(Ms, S)', [], 1);
CompareSecondPartU = reshape(Mu*khatrirao(V, W)' + V*khatrirao(W, Mu)' + W*khatrirao(Mu, V)', [], 1);
CompareSecondPartV = reshape(Mv*khatrirao(W, U)' + W*khatrirao(U, Mv)' + U*khatrirao(Mv, W)', [], 1);
CompareSecondPartW = reshape(Mw*khatrirao(U, V)' + U*khatrirao(V, Mw)' + V*khatrirao(Mw, U)', [], 1);

[norm(Js*Ms(:) - CompareSecondPartS); norm(Ju*Mu(:) - CompareSecondPartU); norm(Jv*Mv(:) - CompareSecondPartV); norm(Jw*Mw(:) - CompareSecondPartW)]

CompareFirstPartS = kron(khatrirao(S, S)', I) + kron(khatrirao(S, S)', I)*P2 + kron(khatrirao(S, S)', I)*P3;
CompareFirstPartU = kron(khatrirao(V, W)', I) + kron(khatrirao(W, V)', I)*P2 + kron(khatrirao(V, W)', I)*P3;
CompareFirstPartV = kron(khatrirao(W, U)', I) + kron(khatrirao(U, W)', I)*P2 + kron(khatrirao(W, U)', I)*P3;
CompareFirstPartW = kron(khatrirao(U, V)', I) + kron(khatrirao(V, U)', I)*P2 + kron(khatrirao(U, V)', I)*P3;

[norm(Js' - CompareFirstPartS); norm(Ju' - CompareFirstPartU); norm(Jv' - CompareFirstPartV); norm(Jw' - CompareFirstPartW)]

Math_TwoB = reshape(CompareFirstPartU*CompareSecondPartU, [], 2)