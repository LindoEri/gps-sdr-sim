
using System.Net.Sockets;

namespace GPS模拟
{
    static class TCP
    {
        static TcpClient t;

        public static bool TCPConnect(string ip, int port)
        {
            try
            {
                t = new TcpClient();
                t.Connect(ip, port);
                return true;
            }
            catch
            {
                return false;
            }
        }
        public static void TCPclose()
        {
            if (t==null || t.Connected == false) return;
            t.Close();
        }
        public static void TCPsend(byte[] datap)
        {
            if (t.Connected == false) return;
            try
            {
                NetworkStream ns = t.GetStream();
                ns.Write(datap, 0, datap.Length);
            }
            catch
            {

            }
        }

    }
}
