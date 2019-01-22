using GMap.NET;
using GMap.NET.MapProviders;
using GMap.NET.WindowsPresentation;
using Microsoft.Win32;
using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.IO;
using System.Net;
using System.Threading;
using System.Windows;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Shapes;

namespace GPS模拟
{
    /// <summary>
    /// MainWindow.xaml 的交互逻辑
    /// </summary>
    public partial class MainWindow : Window
    {
        const double R2D = 57.2957795131;

        public MainWindow()
        {
            InitializeComponent();
        }

        class LocInf : INotifyPropertyChanged
        {
            private double[] Loc;

            public event PropertyChangedEventHandler PropertyChanged;
            public void OnPropertyChanged()
            {
                PropertyChangedEventHandler handler = this.PropertyChanged;
                handler(this, new PropertyChangedEventArgs("Lng"));
                handler(this, new PropertyChangedEventArgs("Lat"));
                handler(this, new PropertyChangedEventArgs("Ele"));
            }

            public string Lng
            {
                get
                {
                    return (Loc[0] * R2D).ToString("0.0000000");
                }
            }
            public string Lat
            {
                get
                {
                    return (Loc[1] * R2D).ToString("0.0000000");
                }
            }
            public string Ele
            {
                get
                {
                    return Loc[2].ToString("0.0000000");
                }
            }
            public LocInf(double[] pos)
            {
                Loc = pos;
            }
        }

        Point tmp = new Point();
        LocInf now, tg;
        Thread mov, th;
        double[] target;
        string nav;
        double DefaultHeight;
        string APIKey = null;

        private void Window_Loaded(object sender, RoutedEventArgs e)
        {
            GPS.Loc = new double[3];
            Dictionary<string, object> Settings = ReadSettings();
            GPS.Loc[0] = (double)Settings["Lng"];
            GPS.Loc[1] = (double)Settings["Lat"];
            GPS.Loc[2] = (double)Settings["Ele"];
            IPAddr.Text = (string)Settings["Addr"];
            Port.Text = ((int)Settings["Port"]).ToString();
            EleDefault.Text = ((double)Settings["Height"]).ToString();
            EleAPIKey.Text = ((string)Settings["EleKey"]).ToString();

            g.MapProvider = GMapProviders.GoogleChinaMap;
            g.Position = new GMap.NET.PointLatLng(GPS.Loc[1] * R2D, GPS.Loc[0] * R2D);
            g.MinZoom = 1;
            g.MaxZoom = 19;
            g.Zoom = (double)Settings["Zoom"];
            GMapProvider.Language = LanguageType.ChineseSimplified;
            g.MouseRightButtonDown += new MouseButtonEventHandler(GetPoint);
            g.OnMapZoomChanged += new MapZoomChanged(DrawMarker);
            g.ShowCenter = false;
            g.DragButton = MouseButton.Left;

            now = new LocInf(GPS.Loc);
            Info.DataContext = now;
            target = new double[3];
            for (int i = 0; i < 3; i++)
                target[i] = GPS.Loc[i];
            tg = new LocInf(target);
            Target.DataContext = tg;
            DrawMarker();
        }
        void GetPoint(object sender, MouseButtonEventArgs e)
        {
            tmp = e.GetPosition((IInputElement)g);
        }
        private void Window_Closing(object sender, System.ComponentModel.CancelEventArgs e)
        {
            TCP.TCPclose();
            if (mov != null && mov.IsAlive == true) mov.Abort();
            if (th != null && th.IsAlive == true) th.Abort();
            string tempFile = System.IO.Path.GetTempFileName();
            if (File.Exists(tempFile + "navfile-gps") == true) File.Delete(tempFile + "navfile-gps");
            if (File.Exists(tempFile + "navfile-gps-d") == true) File.Delete(tempFile + "navfile-gps-d");

            try
            {
                FileStream fs = new FileStream(".\\Settings.ini", FileMode.Create);
                StreamWriter sw = new StreamWriter(fs, System.Text.Encoding.ASCII);
                sw.WriteLine("Longitude in Rad=" + GPS.Loc[0].ToString());
                sw.WriteLine("Latitude in Rad=" + GPS.Loc[1].ToString());
                sw.WriteLine("Elevation=" + GPS.Loc[2].ToString());
                sw.WriteLine("IP Address=" + IPAddr.Text);
                sw.WriteLine("Port=" + Port.Text);
                sw.WriteLine("Zoom Level=" + g.Zoom.ToString());
                sw.WriteLine("Default Height=" + EleDefault.Text);
                sw.WriteLine("Google Elevation Service Key=" + EleAPIKey.Text);
                sw.Close();
                fs.Close();
            }
            catch { }
        }
        private void Connnect_Click(object sender, RoutedEventArgs e)
        {
            try
            {
                if (GetRINEXFile() == false)
                {
                    MessageBox.Show("获取星历失败！", "GPS Simulator", MessageBoxButton.OK, MessageBoxImage.Warning);
                    return;
                }

                if (double.TryParse(EleDefault.Text, out DefaultHeight) == false)
                {
                    MessageBox.Show("纬度默认值设置失败！", "GPS Simulator", MessageBoxButton.OK, MessageBoxImage.Warning);
                    return;
                }

                if (EleAPIKey.Text != "")
                    APIKey = EleAPIKey.Text;

                if (TCP.TCPConnect(IPAddr.Text, int.Parse(Port.Text)) == true)
                {
                    th = new Thread(() =>
                    {
                        GPS.GPS_SIM(nav);
                    });
                    th.Start();

                    Option.Visibility = Visibility.Hidden;
                    Operate.Visibility = Visibility.Visible;

                }
                else
                    MessageBox.Show("连接失败！", "GPS Simulator", MessageBoxButton.OK, MessageBoxImage.Warning);
            }
            catch
            {
                MessageBox.Show("连接失败！", "GPS Simulator", MessageBoxButton.OK, MessageBoxImage.Warning);
            }

        }
        private void LzwUncompress(string path, string unzippedPath)
        {
            byte[] buffer = new byte[4096];

            using (Stream inStream = new Ebixio.LZW.LzwInputStream(File.OpenRead(path)))
            using (FileStream outStream = File.Create(unzippedPath))
            {
                int read;
                while ((read = inStream.Read(buffer, 0, buffer.Length)) > 0)
                {
                    outStream.Write(buffer, 0, read);
                }
            }
        }
        private void Getllh(double[] Out)
        {
            PointLatLng point = g.FromLocalToLatLng((int)tmp.X, (int)tmp.Y);
            double Ele = DefaultHeight;
            if (APIKey != null)
            {
                try
                {
                    WebRequest request = WebRequest.Create(string.Format("https://maps.googleapis.com/maps/api/elevation/json?locations={0},{1}&key=" + APIKey, point.Lat, point.Lng));
                    request.Timeout = 1000;
                    WebResponse response = request.GetResponse();
                    string sr = new StreamReader(response.GetResponseStream() ?? new MemoryStream()).ReadToEnd();
                    int pos = sr.IndexOf("elevation") + "elevation'=".Length + 2;
                    string tmps = sr.Substring(pos);
                    pos = tmps.IndexOf(",");
                    tmps = tmps.Substring(0, pos);
                    Ele = double.Parse(tmps);
                }
                catch
                {
                }
            }

            Out[0] = point.Lng / R2D;
            Out[1] = point.Lat / R2D;
            Out[2] = Ele;
        }
        private void DrawMarker()
        {
            g.Markers.Clear();
            GMapMarker m = new GMapMarker(new PointLatLng(GPS.Loc[1] * R2D, GPS.Loc[0] * R2D));
            m.Shape = new Ellipse
            {

                Width = 5 * Math.Sqrt(g.Zoom),
                Height = 5 * Math.Sqrt(g.Zoom),
                Fill = Brushes.Red,
                Stroke = Brushes.Orange,
                StrokeThickness = 1 * Math.Sqrt(g.Zoom),
            };
            m.Offset = new Point(-3 * Math.Sqrt(g.Zoom), -3 * Math.Sqrt(g.Zoom));
            g.Markers.Add(m);
        }
        private void DirectMove_Click(object sender, RoutedEventArgs e)
        {
            if (mov != null && mov.IsAlive == true) mov.Abort();
            Getllh(GPS.Loc);
            for (int i = 0; i < 3; i++)
                target[i] = GPS.Loc[i];
            g.Position = new PointLatLng(target[1] * R2D, target[0] * R2D);
            now.OnPropertyChanged();
            tg.OnPropertyChanged();
            DrawMarker();
        }
        private void LinearMove_Click(object sender, RoutedEventArgs e)
        {
            Getllh(target);
            tg.OnPropertyChanged();

            double[] vecin2D = new double[3];
            //海拔
            vecin2D[2] = target[2] - GPS.Loc[2];
            //纬度
            vecin2D[1] = (target[1] * R2D - GPS.Loc[1] * R2D) * 110.94 * 1000;
            //经度（设纬度差距不大）
            vecin2D[0] = (target[0] * R2D - GPS.Loc[0] * R2D) * 6371 * 1000 * Math.Cos(vecin2D[0]) / 360;
            double distance = Math.Sqrt(vecin2D[0] * vecin2D[0] + vecin2D[1] * vecin2D[1] + vecin2D[2] * vecin2D[2]);

            double[] vecinllh = new double[3];
            for (int i = 0; i < 3; i++)
                vecinllh[i] = target[i] - GPS.Loc[i];

            if (mov != null && mov.IsAlive == true) mov.Abort();
            mov = new Thread(() =>
            {
                while (true)
                {
                    double v = 0, inaccu = -1;
                    Dispatcher.Invoke(new Action(() => { v = Spd.Value; inaccu = SpdInaccu.Value; }));

                    while (v == 0 || inaccu == -1) ;

                    double PercentinT = v / 3.6 * 0.1 / distance;
                    if (inaccu != 0)
                    {
                        Random r = new Random();
                        PercentinT *= 1 + 2 * (r.NextDouble() - 0.5) * inaccu;
                    }
                    double[] dllh = new double[3];
                    for (int i = 0; i < 3; i++)
                        dllh[i] = PercentinT * vecinllh[i];

                    bool complete = true;
                    for (int i = 0; i < 3; i++)
                    {
                        if (Math.Abs(dllh[i]) >= Math.Abs(target[i] - GPS.Loc[i]))
                            GPS.Loc[i] = target[i];
                        else
                            GPS.Loc[i] += dllh[i];

                        if (target[i] != GPS.Loc[i])
                            complete = false;
                    }


                    Dispatcher.Invoke(new Action(() => { now.OnPropertyChanged(); DrawMarker(); }));
                    if (complete == true) return;
                    Thread.Sleep(100);
                }
            });

            mov.Start();



            now.OnPropertyChanged();
        }
        private void GPSInaccu_ValueChanged(object sender, RoutedPropertyChangedEventArgs<double> e)
        {
            GPS.inacc = GPSInaccu.Value;
        }
        private void FileNav_Checked(object sender, RoutedEventArgs e)
        {
            OpenFileDialog o = new OpenFileDialog();
            if (o.ShowDialog() == true)
                nav = o.FileName;
            else
            {
                e.Handled = true;
                FileNav.IsChecked = false;
                AutoNav.IsChecked = true;
            }
        }
        private void Disconnect_Click(object sender, RoutedEventArgs e)
        {

            TCP.TCPclose();
            if (mov != null && mov.IsAlive == true) mov.Abort();
            if (th != null && th.IsAlive == true) th.Abort();

            Option.Visibility = Visibility.Visible;
            Operate.Visibility = Visibility.Hidden;
        }

        private void Label_MouseDoubleClick(object sender, MouseButtonEventArgs e)
        {
            MessageBox.Show(@"
     LzwInputStream.cs：
    
     Copyright (C) 2009 Gabriel Burca
    
     This program is free software; you can redistribute it and/or
     modify it under the terms of the GNU General Public License
     as published by the Free Software Foundation; either version 2
     of the License, or (at your option) any later version.
    
     This program is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU General Public License for more details.
    
     You should have received a copy of the GNU General Public License
     along with this program; if not, write to the Free Software
     Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
    
     Linking this library statically or dynamically with other modules is
     making a combined work based on this library.  Thus, the terms and
     conditions of the GNU General Public License cover the whole
     combination.


    GPS.cs:
        
    MIT License

    Copyright(c) 2015 Takuji Ebinuma

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (#the Software#), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and / or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

            The above copyright notice and this permission notice shall be included in all
            copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED #AS IS#, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
".Replace('#','"'), "开源许可");

        }

        private Dictionary<string, object> ReadSettings()
        {
            Dictionary<string, object> Settings = new Dictionary<string, object>();
            if (File.Exists(".\\Settings.ini") == true)
            {
                try
                {
                    FileStream fs = new FileStream(".\\Settings.ini", FileMode.Open);
                    StreamReader sr = new StreamReader(fs, System.Text.Encoding.ASCII);
                    string str = sr.ReadLine();
                    while (str != null)
                    {
                        string[] tmp = str.Split('=');
                        if (tmp.Length < 2) continue;
                        switch (tmp[0])
                        {
                            case "Longitude in Rad":
                                double lng;
                                if (double.TryParse(tmp[1], out lng) == true) Settings.Add("Lng", lng);
                                break;
                            case "Latitude in Rad":
                                double lat;
                                if (double.TryParse(tmp[1], out lat) == true) Settings.Add("Lat", lat);
                                break;
                            case "Elevation":
                                double ele;
                                if (double.TryParse(tmp[1], out ele) == true) Settings.Add("Ele", ele);
                                break;
                            case "IP Address":
                                IPAddress addr;
                                if (IPAddress.TryParse(tmp[1], out addr) == true) Settings.Add("Addr", tmp[1]);
                                break;
                            case "Port":
                                int port;
                                if (int.TryParse(tmp[1], out port) == true) Settings.Add("Port", port);
                                break;
                            case "Zoom Level":
                                double zoom;
                                if (double.TryParse(tmp[1], out zoom) == true) Settings.Add("Zoom", zoom);
                                break;
                            case "Default Height":
                                double height;
                                if (double.TryParse(tmp[1], out height) == true) Settings.Add("Height", height);
                                break;
                            case "Google Elevation Service Key":
                                if (tmp[1] != null) Settings.Add("EleKey", tmp[1]);
                                break;
                            default:
                                break;
                        }
                        str = sr.ReadLine();
                    }
                    sr.Close();
                    fs.Close();
                }
                catch
                { }
            }

            if (Settings.ContainsKey("Lng") == false) Settings.Add("Lng", 116.3974 / R2D);
            if (Settings.ContainsKey("Lat") == false) Settings.Add("Lat", 39.9088 / R2D);
            if (Settings.ContainsKey("Ele") == false) Settings.Add("Ele", 50.0);
            if (Settings.ContainsKey("Addr") == false) Settings.Add("Addr", "127.0.0.1");
            if (Settings.ContainsKey("Port") == false) Settings.Add("Port", 3737);
            if (Settings.ContainsKey("Zoom") == false) Settings.Add("Zoom", 12.0);
            if (Settings.ContainsKey("Height") == false) Settings.Add("Height", 100.0);
            if (Settings.ContainsKey("EleKey") == false) Settings.Add("EleKey", "");
            return Settings;
        }
        private bool GetRINEXFile()
        {
            try
            {
                if (AutoNav.IsChecked == true)
                {
                    string path = "ftp://igs.gnsswhu.cn/pub/gps/data/daily/" + DateTime.Now.Year.ToString() + "/";
                    string tempFile = System.IO.Path.GetTempFileName();
                    string ZipFile = tempFile + "navfile-gps";
                    string UnZipFile = tempFile + "navfile-gps-d";

                    bool Flag = false;
                    for (int i = DateTime.Now.DayOfYear; i > 0; i--)
                    {
                        if (Flag == true) break;

                        FtpWebRequest FTP = (FtpWebRequest)FtpWebRequest.Create(new Uri(path + i.ToString()) + "/" + DateTime.Now.Year.ToString().Substring(2) + "n/");
                        FTP.Method = WebRequestMethods.Ftp.ListDirectoryDetails;
                        StreamReader rd = new StreamReader(FTP.GetResponse().GetResponseStream());
                        List<string> fl = new List<string>();
                        while (true)
                        {
                            string s = rd.ReadLine();
                            if (s != null)
                                fl.Add(s);
                            else
                                break;
                        }
                        FTP.Abort();
                        rd.Close();

                        int max = 0;
                        string maxs = "";
                        foreach (string str in fl)
                            if (str.Contains("n.Z") == true)
                            {
                                string strx = str;
                                while (strx.Contains("  ") == true)
                                    strx = strx.Replace("  ", " ");
                                string[] para = strx.Split(' ');
                                if (max < int.Parse(para[4]))
                                {
                                    max = int.Parse(para[4]);
                                    maxs = para[8];
                                }

                            }

                        if (max == 0) continue;

                        FtpWebRequest down = (FtpWebRequest)FtpWebRequest.Create(new Uri(path + i.ToString()) + "/" + DateTime.Now.Year.ToString().Substring(2) + "n/" + maxs);
                        down.UseBinary = true;
                        down.Method = WebRequestMethods.Ftp.DownloadFile;
                        FileStream fs = new FileStream(ZipFile, FileMode.Create);
                        FtpWebResponse fw = (FtpWebResponse)down.GetResponse();
                        Stream sr = fw.GetResponseStream();
                        byte[] buffer = new byte[2048];
                        int bytesRead = sr.Read(buffer, 0, buffer.Length);
                        while (bytesRead > 0)
                        {
                            fs.Write(buffer, 0, bytesRead);
                            bytesRead = sr.Read(buffer, 0, buffer.Length);
                        }
                        fs.Close();
                        Flag = true;
                        break;

                    }
                    LzwUncompress(ZipFile, UnZipFile);
                    nav = UnZipFile;
                }
                return true;
            }
            catch
            { return false; }
        }
    }
}
