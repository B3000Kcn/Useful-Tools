#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# --- Metadata ---
__author__      = "B3000Kcn"
__credits__     = ["DBL1F7E5"]
__copyright__   = "Copyright 2025, B3000Kcn"
__license__     = "MIT"
__version__     = "0.1.0"

"""
smart_layout_brick_labels.py

为 Brick 导出的环形基因组可视化 JSON 执行标签防遮挡优化：
通过「分层错落布局」(Staggered Layout) 自动调整 label 环的 lineLength，
缓解水平标签在圆周上彼此遮挡的问题。

用法示例
--------
1) 最简单用法（输出到自动生成的 *_smart_layout.json）::

    python smart_layout_brick_labels.py -i data_examples/brick_rotated.json

2) 指定输出文件路径::

    python smart_layout_brick_labels.py \\
        -i data_examples/brick_rotated.json \\
        -o data_examples/brick_smart_layout.json

3) 调整基础连线长度、步长与层数（更“张牙舞爪”一点）::

    python smart_layout_brick_labels.py \\
        -i data_examples/brick_rotated.json \\
        -o data_examples/brick_smart_layout.json \\
        -L 60 -S 30 -N 4

参数说明
--------
-i, --input        输入 Brick JSON 文件路径（必选）
-o, --output       输出 JSON 文件路径（可选，不填则在文件名后加 _smart_layout）
-L, --base-length  最内层标签的基础 lineLength，默认 40
-S, --layer-step   每增加一层，lineLength 增加的步长，默认 25
-N, --max-layers   最多使用多少层来错落标签，默认 3

注意事项
--------
- 仅处理 type == "label" 的环；
- 只调整/写入 label.data[*].lineLength 字段，不改变 lineAngle；
- 如标签原本已有 lineLength，则取 max(原值, 计算值)，避免把人工拉长的线缩短。
"""

import json
import os
import argparse
import math
from typing import List, Dict, Any, Tuple


def compute_base_angle(start: float, end: float, total_length: float) -> float:
    """
    根据基因组区间计算标签的基础角度 (0-360, 12 点方向为 0, 顺时针)。
    """
    mid_point = (start + end) / 2.0
    angle = (mid_point / total_length) * 360.0
    # 归一化到 [0, 360)
    angle = angle % 360.0
    return angle


def angle_distance(a: float, b: float) -> float:
    """
    圆周角度差，返回最小差值 (0-180)。
    """
    diff = abs(a - b) % 360.0
    if diff > 180.0:
        diff = 360.0 - diff
    return diff


def estimate_label_angle_width(text: str, base_char_deg: float = 1.2, min_width: float = 4.0) -> float:
    """
    粗略估计一个水平文字标签在圆周方向上占用的角宽度。
    - base_char_deg: 每个字符占用的角度（经验值，可根据实际效果微调）
    - min_width: 最小角宽度，避免过短导致算法过于激进
    """
    if not text:
        return min_width
    width = max(min_width, len(text) * base_char_deg)
    return width


def assign_layers(
    labels: List[Dict[str, Any]],
    max_layers: int = 3,
    base_length: float = 40.0,
    layer_step: float = 25.0,
    safety_factor: float = 1.05,
) -> None:
    """
    核心“分层错落布局”算法：
    - 输入：已经附带 base_angle 和 angle_width 的 label 列表（一个 ring 内的）
    - 输出：为每个 label 设置合适的 lineLength（内层短，外层长），尽量避免重叠。

    思路（贪心 + 扫描线）：
    1. 按 base_angle 从小到大排序。
    2. 针对每个 label，尝试从第 0 层到第 max_layers-1 层寻找一个不与已有标签重叠的位置。
    3. 每一层维护一个“当前已经占用的角区段终点”（last_end_angle），如果新标签的起点 > last_end_angle，则可以放入本层。
    4. 若所有层都放不下，则放到最外层，并允许一定程度的重叠（保证算法总是可行）。
    """
    # 先按角度排序
    labels.sort(key=lambda x: x["base_angle"])

    # 为每一层维护一个“最后占用到的角度终点”
    # layers[i] = last_end_angle for layer i
    layers: List[float] = [-1e9 for _ in range(max_layers)]

    for label in labels:
        angle = label["base_angle"]
        half_width = (label["angle_width"] * safety_factor) / 2.0
        start_angle = angle - half_width
        end_angle = angle + half_width

        chosen_layer = None

        for layer_index in range(max_layers):
            last_end = layers[layer_index]
            # 考虑圆周连续性问题：我们用一个简单近似——假设标签都比较稀疏，直接用线性比较
            # 如果当前标签起点比该层 last_end 大，就认为不重叠
            if start_angle > last_end:
                chosen_layer = layer_index
                layers[layer_index] = end_angle
                break

        if chosen_layer is None:
            # 所有层都不行，那就硬塞到最外层，并更新最外层的 last_end
            chosen_layer = max_layers - 1
            layers[chosen_layer] = max(layers[chosen_layer], end_angle)

        # 根据层级设置 lineLength：越外层越长
        line_length = base_length + chosen_layer * layer_step
        label["assigned_layer"] = chosen_layer
        label["assigned_lineLength"] = line_length


def smart_layout_labels(
    input_path: str,
    output_path: str = None,
    base_length: float = 40.0,
    layer_step: float = 25.0,
    max_layers: int = 3,
) -> None:
    """
    为 Brick JSON 中的 label 环执行“智能分层错落布局”，缓解标签之间的遮挡问题。

    主要操作：
    1. 读取 JSON，找到 reference 的 total_length。
    2. 找到 type == "label" 的 ring。
    3. 对 ring.data 中每个 label:
       - 计算 base_angle（基于 start/end/total_length）。
       - 粗略估计 angle_width（基于 text 长度）。
    4. 对同一 ring 中的标签运行 assign_layers()，为其分配 layer 和 lineLength。
    5. 将计算结果写回 label["lineLength"]，保留/可选修正 lineAngle。
    6. 保存新的 JSON。
    """
    try:
        with open(input_path, "r", encoding="utf-8") as f:
            data = json.load(f)
    except FileNotFoundError:
        print(f"❌ 错误：找不到文件 {input_path}")
        return

    # 1. 找 reference total_length
    total_length = 0
    for file_obj in data.get("files", []):
        if file_obj.get("type") == "reference":
            total_length = file_obj.get("total_length", 0)
            print(
                f"🔍 找到参考基因组: {file_obj.get('name_original')}，长度: {total_length} bp"
            )
            break

    if total_length == 0:
        print("❌ 错误：未找到参考基因组长度 (total_length)，无法进行智能布局。")
        return

    processed_labels = 0
    ring_count = 0

    # 2. 遍历 rings，找到 label 类型
    for ring in data.get("rings", []):
        if ring.get("type") != "label":
            continue

        ring_count += 1
        ring_title = ring.get("title", f"Label Ring #{ring_count}")
        print(f"\n⚙️ 正在处理标签环: {ring_title} (Index: {ring.get('index')})")

        labels = ring.get("data", [])
        if not labels:
            print("  ⏭️ 该环没有标签数据，跳过。")
            continue

        # 3. 为每个 label 计算 base_angle 和 angle_width
        enriched_labels: List[Dict[str, Any]] = []
        for label in labels:
            start = label.get("start")
            end = label.get("end")
            text = label.get("text", "")

            if start is None or end is None:
                # 没有位置就跳过
                continue

            base_angle = compute_base_angle(float(start), float(end), float(total_length))
            angle_width = estimate_label_angle_width(text)

            enriched = {
                "raw": label,
                "base_angle": base_angle,
                "angle_width": angle_width,
            }
            enriched_labels.append(enriched)

        if not enriched_labels:
            print("  ⚠️ 标签环中没有可用于布局的标签（缺少 start/end）。")
            continue

        # 4. 为本环的标签分配层级和 lineLength
        assign_layers(
            enriched_labels,
            max_layers=max_layers,
            base_length=base_length,
            layer_step=layer_step,
        )

        # 5. 将结果写回原始 label
        for item in enriched_labels:
            label = item["raw"]
            assigned_length = item["assigned_lineLength"]

            # 如果已有 lineLength，取较大值以避免缩短（可按需调整策略）
            existing_length = label.get("lineLength")
            if existing_length is None:
                label["lineLength"] = assigned_length
            else:
                try:
                    # 确保是数值比较
                    existing_val = float(existing_length)
                    label["lineLength"] = max(existing_val, assigned_length)
                except (TypeError, ValueError):
                    # 出现异常就直接覆盖
                    label["lineLength"] = assigned_length

            processed_labels += 1

        print(
            f"  ✅ 本环完成智能布局：共 {len(enriched_labels)} 个标签，使用层数 ≤ {max_layers}"
        )

    # 6. 输出 JSON
    if output_path is None:
        base, ext = os.path.splitext(input_path)
        output_path = f"{base}_smart_layout{ext}"

    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)

    print("\n✅ 全部处理完成！")
    print(f"👀 共调整了 {processed_labels} 个标签的 lineLength。")
    print(f"💾 结果已保存至: {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description=(
            "为 Brick 导出的 JSON 标签执行智能防遮挡布局："
            "基于角度和标签长度自动分层、错落调整 lineLength。"
        )
    )
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="输入 Brick JSON 文件路径",
    )
    parser.add_argument(
        "-o",
        "--output",
        required=False,
        help="输出 JSON 文件路径（可选，不填则在输入文件名后加 _smart_layout）",
    )
    parser.add_argument(
        "-L",
        "--base-length",
        type=float,
        default=40.0,
        help="最内层标签的基础连线长度 (lineLength)，默认 40",
    )
    parser.add_argument(
        "-S",
        "--layer-step",
        type=float,
        default=25.0,
        help="每增加一层，lineLength 增加的长度步长，默认 25",
    )
    parser.add_argument(
        "-N",
        "--max-layers",
        type=int,
        default=3,
        help="最多使用多少层来错落标签，默认 3",
    )

    args = parser.parse_args()

    smart_layout_labels(
        input_path=args.input,
        output_path=args.output,
        base_length=args.base_length,
        layer_step=args.layer_step,
        max_layers=args.max_layers,
    )


if __name__ == "__main__":
    main()
