# 3D特征稳定性及构象保留规则

状态：预登记启动方案，尚未运行。与V2§6和3d_eligibility_rules.json一起使用。
从开发结构按拓扑、MW、RotB、电荷、来源分层选60–100图；保留选择ID和种子。候选50/100/200，每级三个独立seed（17,29,42）；明确记录实际嵌入数、有效优化数、去重数、能窗内数、最终保留数。按V2 ETKDGv3/MMFF94s、RMSD 0.75 Å、相对能窗10 kcal/mol试算。去重按能量升序、并列构象ID排序保留第一代表，重原子对称性匹配RMSD；超过50个时主方案用能量升序前50，固定算法并标记截断，结构多样性选代表另作敏感性。输出保留坐标ID和能量，不能把生成频率当平衡概率。
对每一个拟入模的描述符×统计量（中位数、p10/p90、IQR、min/max、IMHB比例等）单独检验：同seed不同采样规模及同规模独立seed的绝对差/相对差、秩相关；组内失败也计入报告。默认试算容差=max(5%×较高级别参考绝对值, absolute_floor)，绝对floor候选：PSASA/total SASA 2 Å²、Rgyr 0.1 Å、relative SASA/NPR/asphericity/IMHB比例0.05（比例量）、IMHB count统计量0.5、体积5 Å³；各量单位、统计量及适用范围在正式pilot配置登记，不能直接套给未定义特征。建议80%通过为启动门槛，分别报告主要拓扑×MW层及柔性/电荷层，样本过少标记证据不足，不以总体通过掩盖困难层。参考比较采用100与200及200的独立seed；选择最小满足与更高级别一致的构象数；不能证明充分采样时只称操作性稳定。给出每特征每层可计算率、通过率、差值分布、代表数和受截断比例。阈值仅是候选工程设置，pilot结束前在不看预测表现情况下签署正式配置。
不稳定或高缺失的主特征按预定规则去掉整个特征列/降为探索性；不因某个特征失败删D2D测试数据。特征列表锁定后D3D按有限且有效特征完整性冻结；主比较不按测试误差选择特征。开发外推无标签结构用相同采样设置/预登记复杂度升级规则并报告失败率；绝不拿外推标签设定采样。若有困难分子一次追加至400候选（资源允许，规则在正式计算前登记），仍失败标明原因；不用不断扩采直到好结果。未加权描述符为主，min/max未通过则不进入主向量；IMHB fraction只是采样构象比例。L2验证不同近似不等于真实溶液构象验证。
交付：pilot_members.csv、descriptor_stability.csv（graph、feature、statistic、seed、n_requested/n_retained、reference、difference、pass、reason）、stability_by_stratum.csv、feature_lock.json、3d_coverage.csv、conformer_qc.csv、输入输出hash。没有冻结feature_lock不进入正式3D因子研究。
