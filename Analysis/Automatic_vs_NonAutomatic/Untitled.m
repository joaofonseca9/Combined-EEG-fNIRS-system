% %% Statistical analysis
%
% %% Theta Band.
%
% % Auto Uncued vs Cued.
% auto_theta = [autouncued_ERD_ERS_theta autocued_ERD_ERS_theta];
% groups = {'Auto Uncued'; 'Auto Cued'};
% % Test the hypothesis that the data is normaly distributed.
% h_autouncued_theta = adtest(autouncued_ERD_ERS_theta);
% h_autocued_theta = adtest(autocued_ERD_ERS_theta);
% % If normally distributed - ANOVA test.
% if h_autouncued_theta==0 && h_autocued_theta==0
%     [p, tbl, stats] = anova1(auto_theta, groups, 'on');
%     figure;
%     multcompare(stats);
% % If not normally distributed - Friedman's test.
% else
%     [p, tbl, stats] = friedman(auto_theta, 2, 'on');
%     figure;
%     multcompare(stats);
% end
%
% % Non-Auto Uncued vs Cued.
% nonauto_theta = [nonautouncued_ERD_ERS_theta nonautocued_ERD_ERS_theta];
% groups = {'Non-Auto Uncued'; 'Non-Auto Cued'};
% % Test the hypothesis that the data is normaly distributed.
% h_nonautouncued_theta = adtest(nonautouncued_ERD_ERS_theta);
% h_nonautocued_theta = adtest(nonautocued_ERD_ERS_theta);
% % If normally distributed - ANOVA test.
% if h_nonautouncued_theta==0 && h_nonautocued_theta==0
%     [p, tbl, stats] = anova1(nonauto_theta, groups, 'on');
%     figure;
%     multcompare(stats);
% % If not normally distributed - Friedman's test.
% else
%     [p, tbl, stats] = friedman(nonauto_theta, 2, 'on');
%     figure;
%     multcompare(stats);
% end
%
% disp('Theta');
% pause;
%
% %% Alpha Band.
%
% % Auto Uncued vs Cued.
% auto_alpha = [autouncued_ERD_ERS_alpha autocued_ERD_ERS_alpha];
% groups = {'Auto Uncued'; 'Auto Cued'};
% % Test the hypothesis that the data is normaly distributed.
% h_autouncued_alpha = adtest(autouncued_ERD_ERS_alpha);
% h_autocued_alpha = adtest(autocued_ERD_ERS_alpha);
% % If normally distributed - ANOVA test.
% if h_autouncued_alpha==0 && h_autocued_alpha==0
%     [p, tbl, stats] = anova1(auto_alpha, groups, 'on');
%     figure;
%     multcompare(stats);
% % If not normally distributed - Friedman's test.
% else
%     [p, tbl, stats] = friedman(auto_alpha, 2, 'on');
%     figure;
%     multcompare(stats);
% end
%
% % Non-Auto Uncued vs Cued.
% nonauto_alpha = [nonautouncued_ERD_ERS_alpha nonautocued_ERD_ERS_alpha];
% groups = {'Non-Auto Uncued'; 'Non-Auto Cued'};
% % Test the hypothesis that the data is normaly distributed.
% h_nonautouncued_alpha = adtest(nonautouncued_ERD_ERS_alpha);
% h_nonautocued_alpha = adtest(nonautocued_ERD_ERS_alpha);
% % If normally distributed - ANOVA test.
% if h_nonautouncued_alpha==0 && h_nonautocued_alpha==0
%     [p, tbl, stats] = anova1(nonauto_alpha, groups, 'on');
%     figure;
%     multcompare(stats);
% % If not normally distributed - Friedman's test.
% else
%     [p, tbl, stats] = friedman(nonauto_alpha, 2, 'on');
%     figure;
%     multcompare(stats);
% end
%
% disp('Alpha');
% pause;
%
% %% Beta Band.
%
% % Auto Uncued vs Cued.
% auto_beta = [autouncued_ERD_ERS_beta autocued_ERD_ERS_beta];
% groups = {'Auto Uncued'; 'Auto Cued'};
% % Test the hypothesis that the data is normaly distributed.
% h_autouncued_beta = adtest(autouncued_ERD_ERS_beta);
% h_autocued_beta = adtest(autocued_ERD_ERS_beta);
% % If normally distributed - ANOVA test.
% if h_autouncued_beta==0 && h_autocued_beta==0
%     [p, tbl, stats] = anova1(auto_beta, groups, 'on');
%     figure;
%     multcompare(stats);
% % If not normally distributed - Friedman's test.
% else
%     [p, tbl, stats] = friedman(auto_beta, 2, 'on');
%     figure;
%     multcompare(stats);
% end
%
% % Non-Auto Uncued vs Cued.
% nonauto_beta = [nonautouncued_ERD_ERS_beta nonautocued_ERD_ERS_beta];
% groups = {'Non-Auto Uncued'; 'Non-Auto Cued'};
% % Test the hypothesis that the data is normaly distributed.
% h_nonautouncued_beta = adtest(nonautouncued_ERD_ERS_beta);
% h_nonautocued_beta = adtest(nonautocued_ERD_ERS_beta);
% % If normally distributed - ANOVA test.
% if h_nonautouncued_beta==0 && h_nonautocued_beta==0
%     [p, tbl, stats] = anova1(nonauto_beta, groups, 'on');
%     figure;
%     multcompare(stats);
% % If not normally distributed - Friedman's test.
% else
%     [p, tbl, stats] = friedman(nonauto_beta, 2, 'on');
%     figure;
%     multcompare(stats);
% end
%
% disp('Beta');
% pause;
%
% %% Gamma Band.
%
% % Auto Uncued vs Cued.
% auto_gamma = [autouncued_ERD_ERS_gamma autocued_ERD_ERS_gamma];
% groups = {'Auto Uncued'; 'Auto Cued'};
% % Test the hypothesis that the data is normaly distributed.
% h_autouncued_gamma = adtest(autouncued_ERD_ERS_gamma);
% h_autocued_gamma = adtest(autocued_ERD_ERS_gamma);
% % If normally distributed - ANOVA test.
% if h_autouncued_gamma==0 && h_autocued_gamma==0
%     [p, tbl, stats] = anova1(auto_gamma, groups, 'on');
%     figure;
%     multcompare(stats);
% % If not normally distributed - Friedman's test.
% else
%     [p, tbl, stats] = friedman(auto_gamma, 2, 'on');
%     figure;
%     multcompare(stats);
% end
%
% % Non-Auto Uncued vs Cued.
% nonauto_gamma = [nonautouncued_ERD_ERS_gamma nonautocued_ERD_ERS_gamma];
% groups = {'Non-Auto Uncued'; 'Non-Auto Cued'};
% % Test the hypothesis that the data is normaly distributed.
% h_nonautouncued_gamma = adtest(nonautouncued_ERD_ERS_gamma);
% h_nonautocued_gamma = adtest(nonautocued_ERD_ERS_gamma);
% % If normally distributed - ANOVA test.
% if h_nonautouncued_gamma==0 && h_nonautocued_gamma==0
%     [p, tbl, stats] = anova1(nonauto_gamma, groups, 'on');
%     figure;
%     multcompare(stats);
% % If not normally distributed - Friedman's test.
% else
%     [p, tbl, stats] = friedman(nonauto_gamma, 2, 'on');
%     figure;
%     multcompare(stats);
% end
%
% disp('Gamma');
% pause;




% %% Statistical analysis.
%
% %% F7.
%
% % Auto Uncued vs Cued.
% auto_F7 = [autouncued_power_F7 autocued_power_F7];
% groups = {'Auto Uncued'; 'Auto Cued'};
% % Test the hypothesis that the data is normaly distributed.
% h_autouncued_F7 = adtest(autouncued_power_F7);
% h_autocued_F7 = adtest(autocued_power_F7);
% % If normally distributed - ANOVA test.
% if h_autouncued_F7==0 && h_autocued_F7==0
%     [p, tbl, stats] = anova1(auto_F7, groups, 'on');
%     figure;
%     multcompare(stats);
% % If not normally distributed - Friedman's test.
% else
%     [p, tbl, stats] = friedman(auto_F7, 2, 'on');
%     figure;
%     multcompare(stats);
% end
%
% % Non-Auto Uncued vs Non-Cued.
% nonauto_F7 = [nonautouncued_power_F7 nonautocued_power_F7];
% groups = {'Auto Uncued'; 'Auto Cued'};
% % Test the hypothesis that the data is normaly distributed.
% h_autouncued_F7 = adtest(nonautouncued_power_F7);
% h_nonautocued_F7 = adtest(nonautocued_power_F7);
% % If normally distributed - ANOVA test.
% if h_nonautouncued_F7==0 && h_nonautocued_F7==0
%     [p, tbl, stats] = anova1(nonauto_F7, groups, 'on');
%     figure;
%     multcompare(stats);
% % If not normally distributed - Friedman's test.
% else
%     [p, tbl, stats] = friedman(nonauto_F7, 2, 'on');
%     figure;
%     multcompare(stats);
% end
%
% disp('F7');
% pause;
%
% %%
%
% % a = mean(autouncued_power_allSubjects, 1, 'omitnan');
% % b = mean(autocued_power_allSubjects, 1, 'omitnan');
% %
% % a_F7 = a(1, F8_loc, :);
% % b_F7 = b(1, F8_loc, :);
% %
% % a_F7_final = [a_F7(1, 1, 1); a_F7(1, 1, 2)];
% % b_F7_final = [b_F7(1, 1, 1); b_F7(1, 1, 2)];
% %
% % % Auto Uncued vs Cued for F7.
% % auto_F7 = [a_F7_final b_F7_final];
% % groups = {'Auto Uncued'; 'Auto Cued'};
% % [p, tbl, stats] = anova1(auto_F7, groups, 'on');