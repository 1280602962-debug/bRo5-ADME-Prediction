# C2：同样本 D-MPNN 联合训练与供体域微调

适用：Caco-2 D2D原SOURCE_CONNECTIVITY外折及主容量预算；最终留出独立另跑。均用同一V2.1二维图+理化向量、图等权标签、加权MSE；不加3D、不跨终点。用原B/S/M/N同fold/seed成员；以graph及连接关系审核集合，B与供体互斥。M/N沿用同MW分箱。所有策略相同测试成员与指标；B→M为供体域微调，不是目标域微调。

| ID | 训练集合/顺序 | 作用 |
| --- | --- | --- |
| TL-B | B从头训练 | 基础及共享预训练权重 |
| TL-N0 | N从头训练 | 目标域小样本基线 |
| TL-JS | B+S随机初始化联合训练 | S联合 |
| TL-JM | B+M随机初始化联合训练 | M联合 |
| TL-JN | B+N随机初始化联合训练 | N联合 |
| TL-FS | B训练→S head-only | S微调 |
| TL-FM | B训练→M head-only | M微调 |
| TL-FN | B训练→N head-only | N微调 |

关键TL01 = MAE(TL-FN)-MAE(TL-JN)。辅助：各供体F-J、F-B、TL-FN-TL-N0；负迁移以F-J或F-B为明确不同基准，正值为相对指定基准误差增大，报告CI并保留负结果。后续源域B+M vs B+S→相同N目标域微调需新增共同N预算与同样本联合控制；本轮不开展。不将头部微调泛称所有迁移算法。

## 实现控制
1. 预先冻结同一网络架构、损失、批大小、特征维度、随机初始化、训练权重及优化器。建议batch32、AdamW；初始lr=1e-3、微调lr=1e-4仅为pilot候选，须无外层标签地锁定。共享B拟合的2D/target scaler用于B、J、F（target均逆变换），N0仅用N拟合，单列；如果原C01联合采用B+X缩放或其他配置，需专用重训，不能冒称直接复用。不得用供体标签决定B阶段scaler。图编码器及其归一化/缓冲区在head-only阶段冻结并置eval；FFN回归头全部可训练，优化器重置，保存参数hash、可训练参数数、训练模式。B权重在同fold/seed/训练数据/config完全一致时共享，J从随机初始化。
2. 主协议采用固定更新步数，而非以外层test早停。pilot前设候选B=50 epoch、供体head=30 epoch（均在V2.1最多100 epoch/阶段内），以实际batch/epoch和训练侧标签行数确定每对总updates U_B+U_X；对应联合训练使用相同总optimizer updates，按图权重批次计划迭代，不能额外利用其他样本。N0用与X微调相同供体updates并记录。累计独立分子、不同标签数、每图抽样曝光次数、可训练参数、梯度更新数、GPU耗时分别记录。匹配总updates并不匹配每图曝光或计算量；头部冻结和阶段性曝光属于策略差异，不归因于单独预训练机制。若要分解，增加replay/全网络微调须新登记且在资源允许下实施。
3. 首轮不做每臂超参网格；采用预登记固定配置/epochs或一次开发侧嵌套pilot，所有臂相同规则和总选择预算；pilot只能在原外层训练侧内选择，不能共享某一外层测试标签信息到另折参数。若内验证缺乏目标域样本（例如B/S/M），使用预登记固定步数，不借额外N标签。先锁定候选配置、步数及允许的GPU上限；外层早停禁止。训练成本差异完整披露。
4. 任务预算：8策略×最多50 fold-seed=400逻辑任务；原C01中B及3个联合模型只有hash完全等同才能复用，理论下限200新增任务。专用B+三联合重训最多200，故首轮新增最大400条外层任务（不含pilot和技术重试）；每任务技术重试最多1次且仅异常，不因性能低重试。Full fine-tuning最多作为预登记N臂敏感性+50任务，其他供体不全展开；本轮不默认执行。
5. 最终封存阶段仅重新训练锁定的B、JN、FN以及按资源登记的辅助策略，用完整开发池仍按锁定预算抽取；同测试集、同缩放及更新规则，另计资源。

交付：tl_config_lock.json、member_pair_check.csv、scaler_manifest.csv、checkpoint_manifest.csv、training_step_exposure.csv（fold/seed/arm/graph_id/label exposure/update count）、tl_oof.csv.gz、tl_paired_deltas.csv、resource_usage.csv、negative_transfer_summary.csv；失败状态及重试记录。实际配置、数据和计算未完成。
