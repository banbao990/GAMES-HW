using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Rigid_Bunny_by_Shape_Matching : MonoBehaviour {
    public bool Launched = false;
    Vector3[] X;
    Vector3[] Q;                            // ri
    Vector3[] V;                            // vi(速度)
    Matrix4x4 QQt = Matrix4x4.zero;         // \sum{r_i * r_i^t}


    float Dt = 0.015f;
    Vector3 g = new Vector3(0, -9.8f, 0);   // 重力加速度
    float mu_n = 0.5f;
    float mu_t = 0.9f;

    // Start is called before the first frame update
    void Start() {
        Mesh mesh = GetComponent<MeshFilter>().mesh;
        V = new Vector3[mesh.vertices.Length];
        X = mesh.vertices;
        Q = mesh.vertices;

        // Centerizing Q.
        Vector3 c = Vector3.zero;
        for (int i = 0; i < Q.Length; i++)
            c += Q[i];
        c /= Q.Length;
        for (int i = 0; i < Q.Length; i++)
            Q[i] -= c;

        // Get QQ^t ready.
        for (int i = 0; i < Q.Length; i++) {
            QQt[0, 0] += Q[i][0] * Q[i][0];
            QQt[0, 1] += Q[i][0] * Q[i][1];
            QQt[0, 2] += Q[i][0] * Q[i][2];
            QQt[1, 0] += Q[i][1] * Q[i][0];
            QQt[1, 1] += Q[i][1] * Q[i][1];
            QQt[1, 2] += Q[i][1] * Q[i][2];
            QQt[2, 0] += Q[i][2] * Q[i][0];
            QQt[2, 1] += Q[i][2] * Q[i][1];
            QQt[2, 2] += Q[i][2] * Q[i][2];
        }
        QQt[3, 3] = 1;

        for (int i = 0; i < X.Length; i++)
            V[i][0] = 4.0f;

        Update_Mesh(transform.position, Matrix4x4.Rotate(transform.rotation), 0);
        transform.position = Vector3.zero;
        transform.rotation = Quaternion.identity;
    }

    // Polar Decomposition that returns the rotation from F.
    Matrix4x4 Get_Rotation(Matrix4x4 F) {
        Matrix4x4 C = Matrix4x4.zero;
        for (int ii = 0; ii < 3; ii++)
            for (int jj = 0; jj < 3; jj++)
                for (int kk = 0; kk < 3; kk++)
                    C[ii, jj] += F[kk, ii] * F[kk, jj];

        Matrix4x4 C2 = Matrix4x4.zero;
        for (int ii = 0; ii < 3; ii++)
            for (int jj = 0; jj < 3; jj++)
                for (int kk = 0; kk < 3; kk++)
                    C2[ii, jj] += C[ii, kk] * C[jj, kk];

        float det = F[0, 0] * F[1, 1] * F[2, 2] +
                        F[0, 1] * F[1, 2] * F[2, 0] +
                        F[1, 0] * F[2, 1] * F[0, 2] -
                        F[0, 2] * F[1, 1] * F[2, 0] -
                        F[0, 1] * F[1, 0] * F[2, 2] -
                        F[0, 0] * F[1, 2] * F[2, 1];

        float I_c = C[0, 0] + C[1, 1] + C[2, 2];
        float I_c2 = I_c * I_c;
        float II_c = 0.5f * (I_c2 - C2[0, 0] - C2[1, 1] - C2[2, 2]);
        float III_c = det * det;
        float k = I_c2 - 3 * II_c;

        Matrix4x4 inv_U = Matrix4x4.zero;
        if (k < 1e-10f) {
            float inv_lambda = 1 / Mathf.Sqrt(I_c / 3);
            inv_U[0, 0] = inv_lambda;
            inv_U[1, 1] = inv_lambda;
            inv_U[2, 2] = inv_lambda;
        } else {
            float l = I_c * (I_c * I_c - 4.5f * II_c) + 13.5f * III_c;
            float k_root = Mathf.Sqrt(k);
            float value = l / (k * k_root);
            if (value < -1.0f) value = -1.0f;
            if (value > 1.0f) value = 1.0f;
            float phi = Mathf.Acos(value);
            float lambda2 = (I_c + 2 * k_root * Mathf.Cos(phi / 3)) / 3.0f;
            float lambda = Mathf.Sqrt(lambda2);

            float III_u = Mathf.Sqrt(III_c);
            if (det < 0) III_u = -III_u;
            float I_u = lambda + Mathf.Sqrt(-lambda2 + I_c + 2 * III_u / lambda);
            float II_u = (I_u * I_u - I_c) * 0.5f;


            float inv_rate, factor;
            inv_rate = 1 / (I_u * II_u - III_u);
            factor = I_u * III_u * inv_rate;

            Matrix4x4 U = Matrix4x4.zero;
            U[0, 0] = factor;
            U[1, 1] = factor;
            U[2, 2] = factor;

            factor = (I_u * I_u - II_u) * inv_rate;
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    U[i, j] += factor * C[i, j] - inv_rate * C2[i, j];

            inv_rate = 1 / III_u;
            factor = II_u * inv_rate;
            inv_U[0, 0] = factor;
            inv_U[1, 1] = factor;
            inv_U[2, 2] = factor;

            factor = -I_u * inv_rate;
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    inv_U[i, j] += factor * U[i, j] + inv_rate * C[i, j];
        }

        Matrix4x4 R = Matrix4x4.zero;
        for (int ii = 0; ii < 3; ii++)
            for (int jj = 0; jj < 3; jj++)
                for (int kk = 0; kk < 3; kk++)
                    R[ii, jj] += F[ii, kk] * inv_U[kk, jj];
        R[3, 3] = 1;
        return R;
    }

    // Update the mesh vertices according to translation c and rotation R.
    // It also updates the velocity.
    void Update_Mesh(Vector3 c, Matrix4x4 R, float inv_dt) {
        for (int i = 0; i < Q.Length; i++) {
            Vector3 x = (Vector3)(R * Q[i]) + c;

            V[i] += (x - X[i]) * inv_dt;
            X[i] = x;
        }
        Mesh mesh = GetComponent<MeshFilter>().mesh;
        mesh.vertices = X;
    }

    void Collision(float inv_dt) {
        Collision_Impulse(new Vector3(0, 0.01f, 0), new Vector3(0, 1, 0));
        Collision_Impulse(new Vector3(2, 0, 0), new Vector3(-1, 0, 0));
    }

    void Collision_Impulse(Vector3 P, Vector3 N) {
        // 直接对每一个点做更新
        int length = Q.Length;
        for (int i = 0; i < length; ++i) {

            // (1) 更新位置
            float phi = Vector3.Dot(X[i] - P, N);
            // 位置在平面上方
            if (phi >= 0) {
                continue;
            }
            X[i] -= phi * N;

            // (2) 更新速度
            float vDotN = Vector3.Dot(V[i], N);
            // 速度和法向同向
            if (vDotN >= 0) {
                continue;
            }

            // 更新
            Vector3 vn = vDotN * N;
            Vector3 vt = V[i] - vn;
            float a = Mathf.Max(0, 1 - mu_t * (1 + mu_n) * vn.magnitude / vt.magnitude);
            vn = -mu_n * vn;
            vt = a * vt;
            V[i] = vt + vn;
        }
    }

    // Update is called once per frame
    void Update() {
        // Game Control
        if (Input.GetKey("r")) {
            Launched = false;
        }
        if (Input.GetKey("l")) {
            Launched = true;
        }
        if (!Launched) { return; }


        //////////////////////////////////////////////
        // Step 1: run a simple particle system.
        //////////////////////////////////////////////

        // 简单更新每一个点的位置
        for (int i = 0; i < V.Length; i++) {
            // 只有重力
            V[i] += Dt * g;
            X[i] += Dt * V[i];
        }

        //////////////////////////////////////////////
        // Step 2: Perform simple particle collision.
        //////////////////////////////////////////////

        Collision(1 / Dt);

        //////////////////////////////////////////////
        // Step 3: Use shape matching to get new translation c and 
        // new rotation R. Update the mesh by c and R.
        //////////////////////////////////////////////

        // Shape Matching (translation)
        Vector3 c = Vector3.zero;
        int length = Q.Length;
        foreach (Vector3 i in X) {
            c += i;
        }
        c /= length;

        // Shape Matching (rotation)
        Matrix4x4 A = Matrix4x4.zero;
        for (int i = 0; i < length; ++i) {
            A = a_plus_b(a_times_b_transpose(X[i] - c, Q[i]), A);
        }
        A[3, 3] = 1;
        A = A * QQt.inverse;
        Matrix4x4 R = Get_Rotation(A);

        Update_Mesh(c, R, 1 / Dt);
    }



    // Utils
    Matrix4x4 Get_Cross_Matrix(Vector3 a) {
        //Get the cross product matrix of vector a
        Matrix4x4 A = Matrix4x4.zero;
        A[0, 0] = 0;
        A[0, 1] = -a[2];
        A[0, 2] = a[1];
        A[1, 0] = a[2];
        A[1, 1] = 0;
        A[1, 2] = -a[0];
        A[2, 0] = -a[1];
        A[2, 1] = a[0];
        A[2, 2] = 0;
        A[3, 3] = 1;
        return A;
    }

    Matrix4x4 a_plus_b(Matrix4x4 a, Matrix4x4 b) {
        Matrix4x4 A = Matrix4x4.zero;
        for (int i = 0; i < 16; ++i) {
            A[i] = a[i] + b[i];
        }
        return A;
    }

    Matrix4x4 a_times_b_transpose(Vector3 a, Vector3 b) {
        Matrix4x4 A = Matrix4x4.zero;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                A[i, j] = a[i] * b[j];
            }
        }
        A[3, 3] = 1;
        return A;
    }
}
