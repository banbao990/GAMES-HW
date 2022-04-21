using UnityEngine;
using System.Collections;
using System;
using System.Collections.Generic;

public class wave_motion : MonoBehaviour {
    const int size = 100;
    const float rate = 0.005f;
    const float rate_inv = 1.0F / rate;
    const float gamma = 0.004f;
    const float gamma_rate = gamma * rate;
    const float damping = 0.996f;
    float[,] old_h;
    float[,] low_h;
    float[,] vh;
    float[,] b;

    bool[,] cg_mask;
    float[,] cg_p;
    float[,] cg_r;
    float[,] cg_Ap;
    // bool tag = true;

    Vector3[] cube_v;
    Vector3[] cube_w;

    const float mass = 10F;
    const float dt = 0.004F;
    const float rho_g_S = 1000F; // \rho_{液} g S_{底面积}
    const float max_water_drop_r = 0.5F;
    System.Random rdm = new System.Random();
    Vector2Int[] delta8;
    Vector2Int[] delta4;

    // Use this for initialization
    void Start() {
        Mesh mesh = GetComponent<MeshFilter>().mesh;
        mesh.Clear();

        Vector3[] X = new Vector3[size * size];

        for (int i = 0; i < size; i++)
            for (int j = 0; j < size; j++) {
                X[i * size + j].x = i * 0.1f - size * 0.05f;
                X[i * size + j].y = 0;
                X[i * size + j].z = j * 0.1f - size * 0.05f;
            }

        int[] T = new int[(size - 1) * (size - 1) * 6];
        int index = 0;
        for (int i = 0; i < size - 1; i++)
            for (int j = 0; j < size - 1; j++) {
                T[index * 6 + 0] = (i + 0) * size + (j + 0);
                T[index * 6 + 1] = (i + 0) * size + (j + 1);
                T[index * 6 + 2] = (i + 1) * size + (j + 1);
                T[index * 6 + 3] = (i + 0) * size + (j + 0);
                T[index * 6 + 4] = (i + 1) * size + (j + 1);
                T[index * 6 + 5] = (i + 1) * size + (j + 0);
                index++;
            }
        mesh.vertices = X;
        mesh.triangles = T;
        mesh.RecalculateNormals();

        low_h = new float[size, size];
        old_h = new float[size, size];
        vh = new float[size, size];
        b = new float[size, size];

        cg_mask = new bool[size, size];
        cg_p = new float[size, size];
        cg_r = new float[size, size];
        cg_Ap = new float[size, size];

        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                low_h[i, j] = 99999;
                old_h[i, j] = 0;
                vh[i, j] = 0;
            }
        }

        cube_v = new Vector3[] { Vector3.zero, Vector3.zero, };
        cube_w = new Vector3[] { Vector3.zero, Vector3.zero, };

        // init delta
        {
            delta8 = new Vector2Int[8];
            int cnt = 0;
            for (int i = -1; i <= 1; ++i) {
                for (int j = -1; j <= 1; ++j) {
                    if (i == 0 && j == 0) { continue; }
                    delta8[cnt++] = new Vector2Int(i, j);
                }
            }

            delta4 = new Vector2Int[4];
            cnt = 0;
            for (int i = -1; i <= 1; ++i) {
                for (int j = -1; j <= 1; ++j) {
                    if ((i == 0 && j == 0) || (i != 0 && j != 0)) { continue; }
                    delta4[cnt++] = new Vector2Int(i, j);
                }
            }
        }
    }

    void A_Times(bool[,] mask, float[,] x, float[,] Ax, int li, int ui, int lj, int uj) {
        for (int i = li; i <= ui; i++)
            for (int j = lj; j <= uj; j++)
                if (i >= 0 && j >= 0 && i < size && j < size && mask[i, j]) {
                    Ax[i, j] = 0;
                    if (i != 0) Ax[i, j] -= x[i - 1, j] - x[i, j];
                    if (i != size - 1) Ax[i, j] -= x[i + 1, j] - x[i, j];
                    if (j != 0) Ax[i, j] -= x[i, j - 1] - x[i, j];
                    if (j != size - 1) Ax[i, j] -= x[i, j + 1] - x[i, j];
                }
    }

    float Dot(bool[,] mask, float[,] x, float[,] y, int li, int ui, int lj, int uj) {
        float ret = 0;
        for (int i = li; i <= ui; i++)
            for (int j = lj; j <= uj; j++)
                if (i >= 0 && j >= 0 && i < size && j < size && mask[i, j]) {
                    ret += x[i, j] * y[i, j];
                }
        return ret;
    }

    void Conjugate_Gradient(bool[,] mask, float[,] b, float[,] x, int li, int ui, int lj, int uj) {
        // Solve the Laplacian problem by CG.
        A_Times(mask, x, cg_r, li, ui, lj, uj);

        for (int i = li; i <= ui; i++)
            for (int j = lj; j <= uj; j++)
                if (i >= 0 && j >= 0 && i < size && j < size && mask[i, j]) {
                    cg_p[i, j] = cg_r[i, j] = b[i, j] - cg_r[i, j];
                }

        float rk_norm = Dot(mask, cg_r, cg_r, li, ui, lj, uj);

        for (int k = 0; k < 128; k++) {
            if (rk_norm < 1e-10f) break;
            A_Times(mask, cg_p, cg_Ap, li, ui, lj, uj);
            float alpha = rk_norm / Dot(mask, cg_p, cg_Ap, li, ui, lj, uj);

            for (int i = li; i <= ui; i++)
                for (int j = lj; j <= uj; j++)
                    if (i >= 0 && j >= 0 && i < size && j < size && mask[i, j]) {
                        x[i, j] += alpha * cg_p[i, j];
                        cg_r[i, j] -= alpha * cg_Ap[i, j];
                    }

            float _rk_norm = Dot(mask, cg_r, cg_r, li, ui, lj, uj);
            float beta = _rk_norm / rk_norm;
            rk_norm = _rk_norm;

            for (int i = li; i <= ui; i++)
                for (int j = lj; j <= uj; j++)
                    if (i >= 0 && j >= 0 && i < size && j < size && mask[i, j]) {
                        cg_p[i, j] = cg_r[i, j] + beta * cg_p[i, j];
                    }
        }

    }

    void Shallow_Wave(float[,] old_h, float[,] h, float[,] new_h) {
        // Step 1:
        // [1.c]: Compute new_h based on the shallow wave model.
        // neumann boundary: 边界高度相等
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                float hij = h[i, j];
                float hij_new = 0;
                int cnt = 0;
                // 向周围 delta[4/8] 个点扩散
                foreach (Vector2Int v in delta4) {
                    int ni = i + v.x, nj = j + v.y;
                    if (ni >= 0 && ni < size && nj >= 0 && nj < size) {
                        ++cnt;
                        hij_new += h[ni, nj];
                    }
                }
                hij_new -= cnt * hij;
                hij_new = hij + (hij - old_h[i, j]) * damping + hij_new * rate;
                new_h[i, j] = hij_new;
            }
        }

        // 每一个 cube 计算完成之后及更新一遍 new_h
        // 和 readme 上不一样, 但是似乎效果好(不需要处理 vh 带来的奇奇怪怪的影响)
        // Step 2: Block->Water coupling
        // [1.d]: for block 1, calculate low_h.
        //      : then set up b and cg_mask for conjugate gradient.
        //      : Solve the Poisson equation to obtain vh (virtual height).
        // [1.d]: for block 2, calculate low_h.
        //      : then set up b and cg_mask for conjugate gradient.
        //      : Solve the Poisson equation to obtain vh (virtual height).
        GameObject[] Cubes = new GameObject[] {
            GameObject.Find("Cube"),
            GameObject.Find("Block"),
        };
        int cube_idx = -1;
        List<Vector3> radius = new List<Vector3>();
        foreach (GameObject cube in Cubes) {
            ++cube_idx;
            radius.Clear();
            // [1.d.1] 获取小方框的包围盒(做了近似: 半径=3)

            // 小方块的中心位置
            Vector3 pos = cube.transform.position;
            Quaternion rot = cube.transform.rotation;
            // 局部坐标下的包围盒
            // Bounds bounds = cube.GetComponent<MeshFilter>().mesh.bounds;
            Bounds bounds = new Bounds(Vector3.zero, Vector3.one); // 和 cube 一样大的 bounds

            // [1.d.2] 从上向下看, 落在包围盒的的液体需要重新计算高度

            // 计算需要计算的位置的下标
            // 之前的计算方式如下
            // X[i * size + j].x = i * 0.1f - size * 0.05f;
            // X[i * size + j].y = 0;
            // X[i * size + j].z = j * 0.1f - size * 0.05f;
            const float tmp_s = size * 0.05F;

            // 边界做一点近似其实没有关系, 左端点+1/右端点-1 保证一定有交点
            // 为了保证有交点, 上面的近似其实还不够(和包围盒有焦点不一定和 cube 有交点)
            // 一个更粗的近似: cube 的内接球在 xOz 平面上的投影的内接 Axis-Aligned 矩形
            // 点的数量近似: (100/10*1)/sqrt(2) = 5sqrt(2) = 7
            // 所以就是中心点向两边各扩展 3 个点(也就是参考答案的写法)
            int il = (int)((pos.x + tmp_s) * 10) - 3;
            int ir = (int)((pos.x + tmp_s) * 10) + 3;
            int jl = (int)((pos.z + tmp_s) * 10) - 3;
            int jr = (int)((pos.z + tmp_s) * 10) + 3;
            il = Mathf.Max(il, 0);
            ir = Mathf.Min(ir, size - 1);
            jl = Mathf.Max(jl, 0);
            jr = Mathf.Min(jr, size - 1);
            // [1.d.3] 计算期望的深度值
            // 参考了参考答案
            // 期望小方块的 cube.y = 0, 根据这一点来计算
            for (int i = il; i <= ir; ++i) {
                for (int j = jl; j <= jr; ++j) {
                    Vector3 p1 = new Vector3(i * 0.1f - tmp_s, 0, j * 0.1f - tmp_s);
                    Vector3 p2 = p1;
                    p1.y = -11F; // 设置一个比较大的值, 保证一定在包围盒的外面
                    p2.y = -10F;
                    // 因为只需要距离, 因此有如下的转化
                    // 光线和 \vec{p2p1} 和现在的 cube (变换后的 cude )求交点
                    // 等价于
                    // 将光线逆变换之后和原来 Axis-Aligned 的 cube 求交
                    // 这样转化的原因是我们可以直接调用 bounds 的和光线求交函数
                    // !! 参考答案里面使用包围盒其实不正确, 应该是使用一个和 cude 一样大的 bound
                    p1 = cube.transform.InverseTransformPoint(p1);
                    p2 = cube.transform.InverseTransformPoint(p2);

                    Ray ray = new Ray(p1, p2 - p1);
                    float dist = 0;
                    bounds.IntersectRay(ray, out dist);
                    // 转化为到 y=0 的距离, 然后因为是负值, 加一个负号
                    low_h[i, j] = dist - 11F; // -(11-dist)

                    radius.Add(p1 + (p2 - p1) * dist - pos);
                }
            }

            // [1.d.4] 共轭梯度下降法
            for (int i = 0; i < size; ++i) {
                for (int j = 0; j < size; ++j) {
                    b[i, j] = 0;
                    vh[i, j] = 0;
                    cg_mask[i, j] = false;
                }
            }
            for (int i = il; i <= ir; ++i) {
                for (int j = jl; j <= jr; ++j) {
                    b[i, j] = (new_h[i, j] - low_h[i, j]) * rate_inv;
                    cg_mask[i, j] = true;
                }
            }
            Conjugate_Gradient(cg_mask, b, vh, il, ir, jl, jr);

            // [1.d.5]: Diminish vh.
            for (int i = 0; i < size; ++i) {
                for (int j = 0; j < size; ++j) {
                    if (cg_mask[i, j]) {
                        vh[i, j] *= gamma_rate;
                    }
                }
            }

            // [1.d.6]: Update new_h by vh.
            for (int i = 0; i < size; ++i) {
                for (int j = 0; j < size; ++j) {
                    // 向周围 delta[4/8] 个点扩散
                    foreach (Vector2Int v in delta4) {
                        int ni = i + v.x, nj = j + v.y;
                        if (ni >= 0 && ni < size && nj >= 0 && nj < size) {
                            new_h[i, j] += vh[ni, nj] - vh[i, j];
                        }
                    }
                }
            }

            // Step 4: Water->Block coupling.
            // 这里实现的效果还是挺烂的
            // 水对 cube 的力放到这里计算
            Vector3 force = Vector3.down * 9.8F * mass;
            Vector3 torque = Vector3.zero;
            int r_idx = 0;
            for (int i = il; i <= ir; ++i) {
                for (int j = jl; j <= jr; ++j) {
                    Vector3 f = Vector3.up * rho_g_S * vh[i, j];
                    force += f;
                    torque += Vector3.Cross( radius[r_idx++], f);
                }
            }
            // 平动
            Vector3 velocity = cube_v[cube_idx];
            velocity = velocity * damping + force / mass * dt;
            cube_v[cube_idx] = velocity;
            cube.transform.position = pos + velocity * dt;
            // 转动
            Vector3 omega = cube_w[cube_idx];
            omega = omega * 0.9F*damping + torque / (50000F * mass); // 看草参考答案的, 不是很懂为什么
            cube_w[cube_idx] = omega;
            const float dt2 = 0.5F * dt;
            Quaternion qt = new Quaternion(omega.x * dt2, omega.y * dt2, omega.z * dt2, 0) * rot;
            qt.x += rot.x;
            qt.y += rot.y;
            qt.z += rot.z;
            qt.w += rot.w;
            cube.transform.rotation = qt;
        }

        // Step 3
        // [1.c]: old_h <- h; h <- new_h;
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                old_h[i, j] = h[i, j];
                h[i, j] = new_h[i, j];
            }
        }
    }

    // Update is called once per frame
    void Update() {
        Mesh mesh = GetComponent<MeshFilter>().mesh;
        Vector3[] X = mesh.vertices;
        float[,] new_h = new float[size, size];
        float[,] h = new float[size, size];

        // [1.a]: Load X.y into h.
        for (int i = 0; i < size; ++i) {
            int ii = i * size;
            for (int j = 0; j < size; ++j) {
                h[i, j] = X[ii + j].y;
            }
        }

        if (Input.GetKeyDown("r")) {
            // TODO: Add random water.
            int ri = rdm.Next(size);
            int rj = rdm.Next(size);
            float r = max_water_drop_r*(float)rdm.NextDouble();
            h[ri, rj] += r;
            // 统计周围有几个点
            int rcnt = 0;
            foreach (Vector2Int v in delta8) {
                int ni = ri + v.x, nj = rj + v.y;
                if (ni >= 0 && ni < size && nj >= 0 && nj < size) {
                    ++rcnt;
                }
            }
            // 周围的液面下降
            float r_each = r / rcnt;
            foreach (Vector2Int v in delta8) {
                int ni = ri + v.x, nj = rj + v.y;
                if (ni >= 0 && ni < size && nj >= 0 && nj < size) {
                    h[ni, nj] -= r_each;
                }
            }
        }

        for (int l = 0; l < 8; l++) {
            Shallow_Wave(old_h, h, new_h);
        }

        // [1.a]: Store h back into X.y and recalculate normal.
        for (int i = 0; i < size; ++i) {
            int ii = i * size;
            for (int j = 0; j < size; ++j) {
                X[ii + j].y = h[i, j];
            }
        }
        mesh.vertices = X;
        mesh.RecalculateNormals();
    }
}