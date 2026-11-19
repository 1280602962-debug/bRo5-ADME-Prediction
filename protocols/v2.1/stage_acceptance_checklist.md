# V2.1 阶段验收（初始全部PENDING）

状态词：PASS / FAIL / PENDING / NOT_ESTIMABLE / STOPPED_WITH_REPORT；仅实测证据允许PASS。每项记录操作人、日期、输入/输出SHA256、代码commit、日志及复核者；不以达到某R²为验收。

| Stage | Required evidence | Stop/degrade |
| --- | --- | --- |
| A | 行级修订/排除、原始record_id、17错标处理、单位复算/守恒、来源映射、Data V1.0 hash、初步3D资格清单 | 错标范围无法界定暂停该终点；不等待所有未知协议 |
| B | 全局封存成员、开发池损耗、source/connectivity/FP06零交叉、预设容量顺序、预算锁定、供体成员/实际分布、完整P01同标签OOF及配对CI/不可估计原因 | 按预设降级；不借封存、不等3D |
| C | model配置、内层选择/早停/权重可复算、实际run matrix与资源上限、同成员OOF和失败日志 | 技术重试一次；未完成不静默删折 |
| D | D3D资格、构象/优化/参数QC、每实际入模特征分层稳定性、锁定feature向量、同原折容量、同成员2D/3D及交互统计 | 不足M/N共同容量则H3 NOT_ESTIMABLE，H1不变 |
| E | 预登记唯一标签和匹配病例、pilot资源/参数/收敛报告，正式病例若通过则分子级响应图与QC | 质量或资源不足停止，不能声称平衡机制 |
| F | 各扩展模块完成/停止报告；模型/特征/家族/人群/脚本/预算签名；留出预测与标签核验、独立统计复算、适用域 | 不反向调参，新增方法探索性或新外部来源 |

实际run matrix由模板展开为 endpoint/population/split/fold/budget/arm/algorithm/representation/donor_seed/model_seed/config_hash/train_members_hash/test_members_hash；缓存命中附相同key证据；每任务记录resource_usage、status、failure_reason、attempt_count。锁定文件列出comparison_registry和各配置/hash。所有失败及效应不明确结果应保留。首轮可立即进入A/B，不需先完成MD资源估算；尚无任何阶段已通过。
