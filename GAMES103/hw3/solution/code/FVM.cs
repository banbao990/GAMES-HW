using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System;
using System.IO;
using System.Linq;

public class FVM : MonoBehaviour {
    const float dt = 0.003f;
    const float mass = 1;
    const float stiffness_0 = 20000.0f;
    const float stiffness_1 = 5000.0f;
    const float damp = 0.999f;

    int[] Tet;
    int tet_number;         // The number of tetrahedra

    Vector3[] Force;
    Vector3[] V;
    Vector3[] X;
    int number;             // The number of vertices

    Matrix4x4[] inv_Dm;

    // For Laplacian smoothing.
    Vector3[] V_sum;
    int[] V_num;

    SVD svd = new SVD();

    // 一些辅助量
    const float EPS = 0.0001F;
    Vector3 Gravity = mass * 9.8F * Vector3.down;   // 重力
    float FloorY = -3F + EPS;
    // 摩擦系数
    float mu_n = 0.5f;
    float mu_t = 0.9f;
    // smooth velocity
    float blend_alpha = 0.8F;

    // Start is called before the first frame update
    void Start() {
        // Application.targetFrameRate = 30; // TODO 提交的时候注释掉这行代码
        // FILO IO: Read the house model from files.
        // The model is from Jonathan Schewchuk's Stellar lib.
        {
            string fileContent = File.ReadAllText("Assets/house2.ele");
            string[] Strings = fileContent.Split(new char[] { ' ', '\t', '\r', '\n' }, StringSplitOptions.RemoveEmptyEntries);

            tet_number = int.Parse(Strings[0]);
            Tet = new int[tet_number * 4];

            for (int tet = 0; tet < tet_number; tet++) {
                Tet[tet * 4 + 0] = int.Parse(Strings[tet * 5 + 4]) - 1;
                Tet[tet * 4 + 1] = int.Parse(Strings[tet * 5 + 5]) - 1;
                Tet[tet * 4 + 2] = int.Parse(Strings[tet * 5 + 6]) - 1;
                Tet[tet * 4 + 3] = int.Parse(Strings[tet * 5 + 7]) - 1;
            }
        }
        {
            string fileContent = File.ReadAllText("Assets/house2.node");
            string[] Strings = fileContent.Split(new char[] { ' ', '\t', '\r', '\n' }, StringSplitOptions.RemoveEmptyEntries);
            number = int.Parse(Strings[0]);
            X = new Vector3[number];
            for (int i = 0; i < number; i++) {
                X[i].x = float.Parse(Strings[i * 5 + 5]) * 0.4f;
                X[i].y = float.Parse(Strings[i * 5 + 6]) * 0.4f;
                X[i].z = float.Parse(Strings[i * 5 + 7]) * 0.4f;
            }
            // Centralize the model.
            Vector3 center = Vector3.zero;
            for (int i = 0; i < number; i++) center += X[i];
            center = center / number;
            for (int i = 0; i < number; i++) {
                X[i] -= center;
                float temp = X[i].y;
                X[i].y = X[i].z;
                X[i].z = temp;
            }
        }
        /*tet_number=1;
        Tet = new int[tet_number*4];
        Tet[0]=0;
        Tet[1]=1;
        Tet[2]=2;
        Tet[3]=3;

        number=4;
        X = new Vector3[number];
        V = new Vector3[number];
        Force = new Vector3[number];
        X[0]= new Vector3(0, 0, 0);
        X[1]= new Vector3(1, 0, 0);
        X[2]= new Vector3(0, 1, 0);
        X[3]= new Vector3(0, 0, 1);*/


        // Create triangle mesh.
        Vector3[] vertices = new Vector3[tet_number * 12];
        int vertex_number = 0;
        for (int tet = 0; tet < tet_number; tet++) {
            vertices[vertex_number++] = X[Tet[tet * 4 + 0]];
            vertices[vertex_number++] = X[Tet[tet * 4 + 2]];
            vertices[vertex_number++] = X[Tet[tet * 4 + 1]];

            vertices[vertex_number++] = X[Tet[tet * 4 + 0]];
            vertices[vertex_number++] = X[Tet[tet * 4 + 3]];
            vertices[vertex_number++] = X[Tet[tet * 4 + 2]];

            vertices[vertex_number++] = X[Tet[tet * 4 + 0]];
            vertices[vertex_number++] = X[Tet[tet * 4 + 1]];
            vertices[vertex_number++] = X[Tet[tet * 4 + 3]];

            vertices[vertex_number++] = X[Tet[tet * 4 + 1]];
            vertices[vertex_number++] = X[Tet[tet * 4 + 2]];
            vertices[vertex_number++] = X[Tet[tet * 4 + 3]];
        }

        int[] triangles = new int[tet_number * 12];
        for (int t = 0; t < tet_number * 4; t++) {
            triangles[t * 3 + 0] = t * 3 + 0;
            triangles[t * 3 + 1] = t * 3 + 1;
            triangles[t * 3 + 2] = t * 3 + 2;
        }
        Mesh mesh = GetComponent<MeshFilter>().mesh;
        mesh.vertices = vertices;
        mesh.triangles = triangles;
        mesh.RecalculateNormals();


        V = new Vector3[number];
        Force = new Vector3[number];
        V_sum = new Vector3[number];
        V_num = new int[number];

        // TODO(F): Need to allocate and assign inv_Dm
        inv_Dm = new Matrix4x4[tet_number];
        for (int i = 0; i < tet_number; ++i) {
            inv_Dm[i] = Build_Edge_Matrix(i);
        }
    }

    void Smooth() {
        for (int i = 0; i < number; i++) {
            V_sum[i] = Vector3.zero;
            V_num[i] = 0;
        }

        // 一个 trick, 直接把当前点的速度加上也没问题, 逻辑更简洁, 运行可能更加快

        for (int tet = 0; tet < tet_number; tet++) {
            int idx = tet * 4;
            for (int i = 0; i < 4; ++i) {
                V_num[Tet[idx + i]] += 3;
                Vector3 v = V[Tet[idx + i]];
                for (int j = 0; j < 4; ++j) {
                    if (i == j) { continue; }
                    V_sum[Tet[idx + j]] += v;
                }
            }
        }

        float blend_beta = 1 - blend_alpha;
        for (int i = 0; i < number; i++) {
            V[i] = blend_alpha * V[i] + blend_beta * V_sum[i] / V_num[i]; // 肯定 V_num[i] > 0
            V_num[i] = 0;
        }
    }

    void Set_Matrix4x4_3x3(ref Matrix4x4 mat) {
        for (int i = 0; i < 3; ++i) {
            mat[i, 3] = mat[3, i] = 0;
        }
        mat[3, 3] = 1F;
    }

    Matrix4x4 Build_Edge_Matrix(int tet) {
        // TODO(F): Need to build edge matrix here.
        Matrix4x4 ret = Build_Edge_Matrix_Xx(tet).inverse;
        Set_Matrix4x4_3x3(ref ret);
        return ret;
        // a00 a01 a02 0
        // a10 a11 a12 0
        // a20 a21 a22 0
        //  0   0   0  1
        // 
        // A 表示上面左上角的 3*3 矩阵
        // 
        // A 0            A^-1 0
        // 0 1 的逆矩阵为  0    1
    }

    Matrix4x4 Build_Edge_Matrix_Xx(int tet) {
        Matrix4x4 ret = Matrix4x4.zero;
        int idx = tet * 4;
        Vector3[] x = new Vector3[3];
        for (int i = 0; i < 3; ++i) {
            x[i] = X[Tet[idx]] - X[Tet[idx + i + 1]];
        }
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                ret[j, i] = x[i][j];
            }
        }
        ret[3, 3] = 1F; // important!
        return ret;
    }

    void _Update() {
        // Jump up.
        if (Input.GetKeyDown(KeyCode.Space)) {
            for (int i = 0; i < number; i++)
                V[i].y += 0.2f;
        }

        for (int i = 0; i < number; i++) {
            // TODO(F): Add gravity to Force.
            Force[i] = Gravity;
        }

        for (int tet = 0; tet < tet_number; tet++) {
            // TODO(F): Deformation Gradient
            Matrix4x4 x = Build_Edge_Matrix_Xx(tet);
            Matrix4x4 F = x * inv_Dm[tet];
            // 矩阵 F 还是如下形式
            // X 0 
            // 0 1
            // 注意接下来只用到了左上角 3*3 的矩阵部分, 因此不需要管他

            // TODO(F): Green Strain
            float tr_G = 0;
            Matrix4x4 G = F.transpose * F;
            for (int i = 0; i < 3; ++i) {
                G[i, i] -= 1F;
                for (int j = 0; j < 3; ++j) {
                    G[i, j] *= 0.5F;
                }
                tr_G += G[i, i];
            }

            // TODO(F): Second PK Stress

            Matrix4x4 S = Matrix4x4.zero, P = Matrix4x4.zero;
#if true
            tr_G *= stiffness_0;
            for (int i = 0; i < 3; ++i) {
                S[i, i] = tr_G;
                for (int j = 0; j < 3; ++j) {
                    S[i, j] += 2F * stiffness_1 * G[i, j];
                }
            }
            P = F * S;
#else
            Matrix4x4 svd_U = Matrix4x4.zero, svd_S = Matrix4x4.zero, svd_V = Matrix4x4.zero;
            svd.svd(F, ref svd_U, ref svd_S, ref svd_V);
            Matrix4x4 Diag = Matrix4x4.zero;
            float l0 = svd_S[0, 0],
                  l1 = svd_S[1, 1],
                  l2 = svd_S[2, 2];
            float I    = l0 * l0 + l1 * l1 + l2 * l2;
            float II   = l0 * l0 * l0 * l0 + l1 * l1 * l1 * l1 + l2 * l2 * l2 * l2;
            float III  = l0 * l1 + l1 * l2 + l2 * l0;
            // PPT 上的 stvk 模型不正确(前面的系数应该是 1/8)
            float dI   = 0.25F * stiffness_0 * (I - 3) - stiffness_1 * 0.5F; 
            float dII  = stiffness_1 * 0.25F;
            // float dIII = 0;

            // 后面反正是 0
            Diag[0, 0] = dI * 2F * l0 + dII * 4F * l0 * l0 * l0;// + dIII * 2F * (l1 * l1 + l2 * l2) * l0;
            Diag[1, 1] = dI * 2F * l1 + dII * 4F * l1 * l1 * l1;// + dIII * 2F * (l2 * l2 + l0 * l0) * l1;
            Diag[2, 2] = dI * 2F * l2 + dII * 4F * l2 * l2 * l2;// + dIII * 2F * (l0 * l0 + l1 * l1) * l2;
            P = svd_U * Diag * svd_V.transpose;
#endif

            // TODO(F): Elastic Force
            float coefficient = -1F / (6F * inv_Dm[tet].determinant);
            Matrix4x4 P2 = P * inv_Dm[tet].transpose;

            int idx = tet * 4;
            Vector3 cumu = Vector3.zero;
            for (int node = 0; node < 3; ++node) {
                Vector3 f = Vector3.zero;
                for (int i = 0; i < 3; ++i) {
                    f[i] = P2[i, node];
                }
                f *= coefficient;
                cumu += f;
                Force[Tet[idx + node + 1]] += f;
            }
            Force[Tet[idx]] -= cumu;
        }

        // 平滑速度, 不加的话会直接炸掉(显式积分的不稳定性)
        Smooth();

        for (int i = 0; i < number; i++) {
            // TODO(F): Update X and V here.
            V[i] += dt * Force[i] / mass;
            V[i] *= damp;
            X[i] += dt * V[i];
            // TODO(F): (Particle) collision with floor.
            {
                // 简单地将平面写死
                // [1] 位置移动到平面外
                if (X[i].y >= FloorY) {
                    continue;
                }
                X[i].y = FloorY;

                // [2] 更新速度, vn 朝外并衰减, vt 衰减
                // 理论上 vt = 0
                if (V[i].y >= 0) {
                    continue;
                }
                Vector3 vn = Vector3.zero;
                vn.y = V[i].y;
                Vector3 vt = V[i] - vn;

                float a = 0;
                if (vt.sqrMagnitude > EPS) {
                    a = Mathf.Max(0, 1 - mu_t * (1 + mu_n) * vn.magnitude / vt.magnitude);
                }
                vt = a * vt;
                vn = -mu_n * vn;
                V[i] = vt + vn;
            }
        }
    }

    // Update is called once per frame
    void Update() {
        for (int l = 0; l < 10; l++)
            _Update();

        // Dump the vertex array for rendering.
        Vector3[] vertices = new Vector3[tet_number * 12];
        int vertex_number = 0;
        for (int tet = 0; tet < tet_number; tet++) {
            vertices[vertex_number++] = X[Tet[tet * 4 + 0]];
            vertices[vertex_number++] = X[Tet[tet * 4 + 2]];
            vertices[vertex_number++] = X[Tet[tet * 4 + 1]];
            vertices[vertex_number++] = X[Tet[tet * 4 + 0]];
            vertices[vertex_number++] = X[Tet[tet * 4 + 3]];
            vertices[vertex_number++] = X[Tet[tet * 4 + 2]];
            vertices[vertex_number++] = X[Tet[tet * 4 + 0]];
            vertices[vertex_number++] = X[Tet[tet * 4 + 1]];
            vertices[vertex_number++] = X[Tet[tet * 4 + 3]];
            vertices[vertex_number++] = X[Tet[tet * 4 + 1]];
            vertices[vertex_number++] = X[Tet[tet * 4 + 2]];
            vertices[vertex_number++] = X[Tet[tet * 4 + 3]];
        }
        Mesh mesh = GetComponent<MeshFilter>().mesh;
        mesh.vertices = vertices;
        mesh.RecalculateNormals();
    }
}
