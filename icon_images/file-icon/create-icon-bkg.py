import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.patheffects as path_effects

def create_file_icon_base_v2(filename="moleditpy_file_icon_v2.png", size_px=1024):
    """
    MoleditPy用のファイルアイコン（ドキュメント型背景）を生成する。
    修正版: 影の形状を紙の形（右上が折れた形）に合わせました。
    """
    dpi = 100
    fig_size_inch = size_px / dpi
    
    # キャンバス準備（背景透明）
    fig, ax = plt.subplots(figsize=(fig_size_inch, fig_size_inch), dpi=dpi)
    # 影が少しずれる分、描画範囲に余裕を持たせる
    ax.set_xlim(-5, 105)
    ax.set_ylim(-5, 105)
    ax.set_aspect('equal')
    ax.axis('off')

    # --- パラメータ設定 ---
    doc_w = 70
    doc_h = 90
    # キャンバスの中央になるように配置
    doc_x = (100 - doc_w) / 2
    doc_y = (100 - doc_h) / 2
    
    fold_size = 20  # 右上の折れ目のサイズ

    # --- 座標定義（紙の本体形状） ---
    # 影の計算にも使うため、先に定義します
    p1 = (doc_x, doc_y)                             # 左下
    p2 = (doc_x + doc_w, doc_y)                     # 右下
    p3 = (doc_x + doc_w, doc_y + doc_h - fold_size) # 右上（折れ目開始点）
    p4 = (doc_x + doc_w - fold_size, doc_y + doc_h) # 右上（折れ目終了点）
    p5 = (doc_x, doc_y + doc_h)                     # 左上
    
    paper_coords = [p1, p2, p3, p4, p5]

    # --- 1. ドロップシャドウ（影） ---
    # 【修正点】紙の形状座標を元に、少し右下にずらした座標を作成
    shadow_off_x = 3
    shadow_off_y = -3
    
    shadow_coords = [
        (p[0] + shadow_off_x, p[1] + shadow_off_y) for p in paper_coords
    ]

    # 多角形で影を描画（joinstyle='round'で角を少し滑らかに）
    shadow_poly = patches.Polygon(shadow_coords, closed=True,
                                  facecolor="black", alpha=0.2, zorder=0,
                                  joinstyle='round')
    
    # 強いぼかし効果を適用してフワッとさせる
    shadow_poly.set_path_effects([
        path_effects.withStroke(linewidth=15, foreground="black", alpha=0.1)
    ])
    ax.add_patch(shadow_poly)

    # --- 2. 紙の本体 ---
    # 先ほど定義した座標で描画
    paper_poly = patches.Polygon(paper_coords, closed=True, 
                                 facecolor="#F9F9F9", edgecolor="#DDDDDD", linewidth=1, zorder=1,
                                 joinstyle='round') # 本体も角を少し滑らかに
    ax.add_patch(paper_poly)

    # --- 3. 紙のテクスチャ（うっすらとした六角形グリッド） ---
    def draw_faint_hexagons():
        hex_size = 8
        h_step = hex_size * np.sqrt(3)
        v_step = hex_size * 1.5
        
        # 描画範囲を少し広めにとってカバー漏れを防ぐ
        for row in range(int((doc_y-10)//v_step), int((doc_y+doc_h+10)//v_step)):
            for col in range(int((doc_x-10)//h_step), int((doc_x+doc_w+10)//h_step)):
                x_pos = col * h_step
                y_pos = row * v_step
                if row % 2 == 1:
                    x_pos += h_step / 2
                
                angles = np.linspace(0, 2*np.pi, 7)
                x_hex = x_pos + hex_size * np.cos(angles)
                y_hex = y_pos + hex_size * np.sin(angles)
                
                poly = patches.Polygon(np.column_stack((x_hex, y_hex)), 
                                      closed=True, edgecolor="#39CCCC", facecolor='none',
                                      linewidth=1, alpha=0.15, zorder=2)
                poly.set_clip_path(paper_poly)
                ax.add_patch(poly)

    draw_faint_hexagons()

    # --- 4. 折れ耳（ドッグイヤー） ---
    fold_coords = [
        p3, # 右端の折れ目
        (doc_x + doc_w - fold_size, doc_y + doc_h - fold_size), # 内側の点
        p4  # 上端の折れ目
    ]
    fold_poly = patches.Polygon(fold_coords, closed=True,
                                facecolor="#EEEEEE", edgecolor="#CCCCCC", linewidth=1, zorder=3,
                                joinstyle='round')
    
    # 折れ耳部分の小さな影
    fold_shadow_coords = [p3, (doc_x + doc_w - fold_size, doc_y + doc_h - fold_size), p4]
    fold_shadow = patches.Polygon(fold_shadow_coords, closed=True, fc="black", alpha=0.1, zorder=2.5, joinstyle='round')
    fold_shadow.set_path_effects([path_effects.withStroke(linewidth=5, foreground="black", alpha=0.1)])
    
    ax.add_patch(fold_shadow)
    ax.add_patch(fold_poly)

    # --- 5. アクセント（下部の帯とテキスト） ---
    bar_h = 10 
    bar_y = doc_y + 12
    bar_rect = patches.Rectangle((doc_x, bar_y), doc_w, bar_h, 
                                 facecolor="#0074D9", alpha=0.9, zorder=2)
    bar_rect.set_clip_path(paper_poly)
    ax.add_patch(bar_rect)
    
    # テキスト設定
    ax.text(doc_x + doc_w/2, bar_y + bar_h/2, "MoleditPy File", 
            ha='center', va='center', 
            fontsize=32,
            color='white', 
            fontweight='bold', 
            fontname='DejaVu Sans',
            zorder=3)

    # 保存（余白調整）
    plt.savefig(filename, bbox_inches='tight', pad_inches=0.1, transparent=True, dpi=dpi)
    plt.close()
    print(f"修正版ファイルアイコン背景を保存しました: {filename}")

if __name__ == "__main__":
    create_file_icon_base_v2()