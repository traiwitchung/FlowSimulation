%%
%refining

block=refine4;

bigblock=zeros(size(block)*2);

bigblock(1:2:end,1:2:end,1:2:end)=block;

bigblock(1:2:end,1:2:end,2:2:end)=block;                                                                                                                                                                                                                                                                                                                                                        

bigblock(1:2:end,2:2:end,1:2:end)=block;                                                                                                                                                                                                                                                                                                                                                        

bigblock(1:2:end,2:2:end,2:2:end)=block;

bigblock(2:2:end,1:2:end,1:2:end)=block;

bigblock(2:2:end,2:2:end,1:2:end)=block;

bigblock(2:2:end,1:2:end,2:2:end)=block;

bigblock(2:2:end,2:2:end,2:2:end)=block;

refine5 = bigblock;

clear bigblock;
clear block;